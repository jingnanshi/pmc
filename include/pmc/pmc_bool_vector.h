#ifndef PMC_BOOL_VECTOR_H_
#define PMC_BOOL_VECTOR_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace pmc {
/// A bare minimum implementation of a boolean vector.
///
/// This class is recommended in place of std::vector<bool> or std::vector<int> for thread-safety.
/// std::vector<bool> could cause a race condition if it is implemented as a dynamic bitset.
/// std::vector<int> is not memory efficient and misleading as an element can hold other than 0 or 1.
class bool_vector {
public:
  bool_vector(std::size_t size = 0UL, bool value = false)
      : data_(size, to_int(value)) {}

  std::size_t size() const noexcept { return data_.size(); }

  bool empty() const noexcept { return data_.empty(); }

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