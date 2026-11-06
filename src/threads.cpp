#ifndef NOTHREADS
#include "GThreads.h"
#endif

#include "ribbon.h"
#include "out.h"
#include "predictions.h"
#include "threads.h"
#include <atomic>
#include <map>
#include <string>

#include "quill/LogMacros.h"
#include "quill/Logger.h"

extern bool RIBBON_DEBUG;
extern quill::Logger *logger;
extern int32_t gene_no;
extern double total_fragments;
extern double total_frag_len;
extern double total_cov_sum;
extern std::map<std::string, double> frag_len_by_cell;
extern std::map<std::string, double> cov_sum_by_cell;

#ifndef NOTHREADS
// Threading: single producer, multiple consumers
// Main thread/program is always loading the producer

extern GMutex data_mutex; // guards the bundle slots that are ready to load
// indices of bundle slots free to load
extern std::vector<int32_t> clear_data_pool;
extern GConditionVar have_bundles; // wakes a thread when a bundle is queued
extern char bundle_work; // bit 0: bundles still being prepared; bit 1: queued

extern GMutex wait_mutex; // controls threads_waiting (idle threads counter)
extern std::atomic<uint8_t> threads_waiting; // idle worker threads
extern GConditionVar have_threads; // wakes the loader when a thread is idle

extern GConditionVar have_clear; // wakes the loader when a slot is free
extern GMutex queue_mutex; // controls bundle_queue and bundles access
extern GFastMutex print_mutex; // for writing the output to file
extern GFastMutex reading_mutex;
extern GFastMutex count_mutex;
#endif

extern bool no_more_bundles;

namespace ribbon {

#ifndef NOTHREADS
  // Waits for a free bundle slot and returns its index.
  int32_t wait_for_data() {
    int32_t idx = -1;

    // Wait until a bundle slot is free and take it.
    data_mutex.lock();
    while (count_of(clear_data_pool) == 0) {
      have_clear.wait(data_mutex);
    }
    idx = clear_data_pool.back();
    clear_data_pool.pop_back();
    data_mutex.unlock();

    return idx;
  }
#endif // NOTHREADS

  // Tells the workers that no more bundles are coming and waits until every
  // worker thread has noticed.
  void wait_for_bundles() {
#ifndef NOTHREADS
    // Record that the BAM file has no more bundles to load.
    reading_mutex.lock();
    no_more_bundles = true;
    reading_mutex.unlock();

    // Tell the workers that no bundles are being prepared anymore.
    queue_mutex.lock();
    bundle_work &= ~(int32_t)0x01; // clear bit 0;
    queue_mutex.unlock();

    // Wake the idle workers until all of them have finished.
    bool are_threads_waiting = true;
    while (are_threads_waiting) {
      are_threads_waiting =
          (threads_waiting.load(std::memory_order_relaxed) > 0);

      if (are_threads_waiting) {
        have_bundles.notify_all();
        current_thread::sleep_for(1);

        are_threads_waiting =
            (threads_waiting.load(std::memory_order_relaxed) > 0);

        current_thread::sleep_for(1);
      }
    }
#else
    no_more_bundles = true;
#endif
  }

  // Infers the transcripts of one bundle, prints them and frees the bundle
  // for reuse.
  void process_bundle(BundleData *bundle) {
    // Log what we are about to process.
    if (RIBBON_DEBUG) {
      LOG_INFO(logger,
               "bundle {}:{}-{} [{} alignments ({} distinct), {} junctions, "
               "{} guides] begins processing...",
               bundle->ref_seq.c_str(), bundle->start, bundle->end,
               bundle->num_reads, count_of(bundle->reads),
               count_of(bundle->junction), count_of(bundle->guides));
    }

    // Assemble the transcripts of this bundle.
    infer_transcripts(bundle);

    // Decide which predictions are kept, then print them; the lock keeps
    // threads from interleaving output.
    if (count_of(bundle->predictions) > 0) {
      finalize_predictions(bundle);
#ifndef NOTHREADS
      GLockGuard<GFastMutex> lock(print_mutex);
#endif
      gene_no = print_results(bundle, gene_no, bundle->ref_seq);
    }

    // Add this bundle's counts to the global totals.
    if (bundle->num_fragments) {
#ifndef NOTHREADS
      GLockGuard<GFastMutex> lock(count_mutex);
#endif
      total_fragments += bundle->num_fragments;
      total_frag_len += bundle->frag_len;
      total_cov_sum += bundle->sum_cov;
      // The same per cell, by barcode (cell numbers are only valid inside one
      // bundle).
      for (const auto& cell : bundle->cell_name) {
        int32_t c = cell.second;
        if (c < count_of(bundle->cell_frag_len))
          frag_len_by_cell[cell.first] += bundle->cell_frag_len[c];
        if (c < count_of(bundle->cell_sum_cov))
          cov_sum_by_cell[cell.first] += bundle->cell_sum_cov[c];
      }
    }

    // Free the bundle data so the slot can be reused.
    bundle->clear();
  }

  // Hands a loaded bundle to the worker threads, or frees it right away if it
  // has no read alignments.
  void push_bundle(BundleData *bundle,
                   std::vector<BundleData *> *bundle_queue) {
    if (count_of(bundle->reads) > 0) {
#ifndef NOTHREADS

      // Add the bundle to the queue for the workers.
      int32_t queue_count = 0;

      queue_mutex.lock();
      bundle_queue->push_back(bundle);
      bundle_work |= 0x02; // set bit 1 to 1
      queue_count = count_of(*bundle_queue);
      queue_mutex.unlock();

      // Wait for an idle worker and wake it.
      wait_mutex.lock();
      while (threads_waiting.load(std::memory_order_relaxed) == 0) {
        have_threads.wait(wait_mutex);
      }
      wait_mutex.unlock();
      have_bundles.notify_one();

      current_thread::yield();

      // Wait until a worker has taken the bundle from the queue.
      queue_mutex.lock();
      while (count_of(*bundle_queue) == queue_count) {
        queue_mutex.unlock();
        have_bundles.notify_one();
        current_thread::yield();
        queue_mutex.lock();
      }
      queue_mutex.unlock();

#else // NOTHREADS = true

      // Without threads we process the bundle right away.
      process_bundle(bundle);
#endif
    }

    // Free a bundle without alignments right away.
    else {
#ifndef NOTHREADS
      data_mutex.lock();
#endif

      bundle->clear();

#ifndef NOTHREADS
      clear_data_pool.push_back(bundle->idx);
      data_mutex.unlock();
#endif
    }
  }

#ifndef NOTHREADS
  // Worker thread: waits for a queued bundle, processes it and gives the
  // bundle slot back to the loader.
  void worker_thread(GThreadData& td) {
    // Get the bundle queue that the loader shares with us.
    std::vector<BundleData *> *bundle_queue =
        (std::vector<BundleData *> *)td.udata;

    // Wait for a ready bundle in the queue
    queue_mutex.lock(); // enter wait-for-notification loop

    while (bundle_work) {
      // Announce that we are idle and wait until a bundle is queued.
      threads_waiting.fetch_add(1, std::memory_order_relaxed);
      queue_mutex.unlock();
      wait_mutex.lock();
      have_threads.notify_one(); // in case main thread is waiting
      wait_mutex.unlock();
      current_thread::yield();
      queue_mutex.lock();
      while (bundle_work && bundle_queue->empty()) {
        // Unlocks queue_mutex and wait until notified
        // When notified, locks queue_mutex and resume
        have_bundles.wait(queue_mutex);
      }

      // We are no longer idle.
      if (threads_waiting.load(std::memory_order_relaxed) > 0)
        threads_waiting.fetch_sub(1, std::memory_order_relaxed);

      // Take a bundle from the queue if one is ready.
      BundleData *ready_bundle = NULL;
      if ((bundle_work & 0x02) != 0) {
        if (!bundle_queue->empty()) {
          ready_bundle = bundle_queue->back();
          bundle_queue->pop_back();
        }

        // Process the bundle outside the queue lock, then give its slot back.
        if (ready_bundle != NULL) {
          if (bundle_queue->empty())
            bundle_work &= ~(int32_t)0x02; // clear bit 1 (queue is empty)

          queue_mutex.unlock();
          process_bundle(ready_bundle);

          // Return the bundle slot to the loader.
          data_mutex.lock();
          clear_data_pool.push_back(ready_bundle->idx);
          data_mutex.unlock();

          have_clear.notify_one(); // inform main thread
          current_thread::yield();

          queue_mutex.lock();
        }
      }
    }
    queue_mutex.unlock();
  }
#endif // NOTHREADS

} // namespace ribbon
