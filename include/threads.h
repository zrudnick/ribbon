#pragma once

namespace ribbon {
#ifndef NOTHREADS
  int32_t wait_for_data();
#endif

  void wait_for_bundles();

  void process_bundle(BundleData *bundle);

  void push_bundle(BundleData *bundle, std::vector<BundleData *> *bundle_queue);

#ifndef NOTHREADS
  void worker_thread(GThreadData& td);
#endif
} // namespace ribbon
