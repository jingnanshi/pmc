#ifndef PMC_BOOL_VECTOR_H_
#define PMC_BOOL_VECTOR_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace pmc {
class bool_vector {
public:
  bool_vector(std::size_t size = 0UL, bool value = false)
      : data_(size, to_int(value)) {}

  void resize(std::size_t size, bool value = false) {
    data_.resize(size, to_int(value));
  }

  class Reference {
  public:
    Reference(std::uint8_t &byte) : byte_(byte) {}

    operator bool() const noexcept { return byte_ != 0; }

    Reference &operator=(bool value) noexcept {
      byte_ = to_int(value);
      return *this;
    }

  private:
    std::uint8_t &byte_;
  };

  bool operator[](std::size_t i) const noexcept { return data_[i] != 0; }

  Reference operator[](std::size_t i) noexcept { return Reference(data_[i]); }

private:
  static constexpr std::uint8_t to_int(bool value) noexcept {
    return value ? 1 : 0;
  }

  std::vector<std::uint8_t> data_;
};
} // namespace pmc

#endif