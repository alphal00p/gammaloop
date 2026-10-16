#include "seed.hpp"
#include "spqr.hpp"
#include <wasi/api.h>

extern "C" __attribute__((
    import_module("typst_env"),
    import_name("wasm_minimal_protocol_write_args_to_buffer"))) void
typst_receive(void *buffer);
extern "C" __attribute__((
    import_module("typst_env"),
    import_name("wasm_minimal_protocol_send_result_to_host"))) void
typst_send(const void *buffer, int size);
extern "C" void __wasm_call_ctors();

// Typst supplies no filesystem or console. These routines are only retained by
// libc++'s fatal-error diagnostics; attempting I/O must trap rather than
// silently pretending it succeeded. Graph operations never call them.
extern "C" __wasi_errno_t __wasi_fd_close(__wasi_fd_t) { __builtin_trap(); }
extern "C" __wasi_errno_t __wasi_fd_seek(__wasi_fd_t, __wasi_filedelta_t,
                                         __wasi_whence_t, __wasi_filesize_t *) {
  __builtin_trap();
}
extern "C" __wasi_errno_t __wasi_fd_write(__wasi_fd_t, const __wasi_ciovec_t *,
                                          size_t, __wasi_size_t *) {
  __builtin_trap();
}

namespace {
int invoke(int size, bool seed) {
  static bool initialized = false;
  if (!initialized) {
    __wasm_call_ctors();
    initialized = true;
  }
  if (size < 0)
    return 1;
  std::string input(size, '\0');
  typst_receive(input.data());
  const auto request = ec::Json::parse(input, nullptr, false);
  if (request.is_discarded()) {
    const std::string error = "Invalid embedding request JSON";
    typst_send(error.data(), error.size());
    return 1;
  }
  const auto result =
      seed ? ec::initialize(request.at("diagram"), request.value("scale", 2.4),
                            request.value("external_sides", true))
           : ec::decompose(request);
  const auto output = result.dump();
  typst_send(output.data(), output.size());
  return 0;
}
} // namespace
extern "C" __attribute__((used, visibility("default"))) int
decompose(int size) {
  return invoke(size, false);
}
extern "C" __attribute__((used, visibility("default"))) int
initialize(int size) {
  return invoke(size, true);
}

// No JavaScript views exist in the Typst host when linear memory grows.
extern "C" void emscripten_notify_memory_growth(int) {}
extern "C" __wasi_errno_t __wasi_fd_read(__wasi_fd_t, const __wasi_iovec_t *,
                                         size_t, __wasi_size_t *) {
  __builtin_trap();
}
// The plugin has no ambient process environment.
extern "C" __wasi_errno_t __wasi_environ_sizes_get(__wasi_size_t *count,
                                                   __wasi_size_t *size) {
  *count = 0;
  *size = 0;
  return 0;
}
extern "C" __wasi_errno_t __wasi_environ_get(uint8_t **, uint8_t *) {
  return 0;
}
