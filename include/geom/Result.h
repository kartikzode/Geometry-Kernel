#pragma once
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace geom {

enum class ErrorCode { InvalidInput, Degenerate, OutOfRange, NotSupported };

struct Error {
    ErrorCode code;
    std::string message;
};

inline Error invalidInput(std::string msg) { return {ErrorCode::InvalidInput, std::move(msg)}; }

// Value-or-error for operations whose failure needs a reason (factories,
// sheet-metal unfolding). Low-level math returns std::optional instead.
// Accessing value() on an error throws std::bad_variant_access.
template <class T>
class Result {
public:
    template <class U, class = std::enable_if_t<!std::is_same_v<std::decay_t<U>, Error> &&
                                                std::is_constructible_v<T, U&&>>>
    Result(U&& value) : v_(std::in_place_index<0>, std::forward<U>(value)) {}
    Result(Error e) : v_(std::in_place_index<1>, std::move(e)) {}

    bool ok() const { return v_.index() == 0; }
    explicit operator bool() const { return ok(); }

    T& value() & { return std::get<0>(v_); }
    const T& value() const& { return std::get<0>(v_); }
    T&& value() && { return std::get<0>(std::move(v_)); }
    T* operator->() { return &value(); }
    const T* operator->() const { return &value(); }
    T& operator*() & { return value(); }
    const T& operator*() const& { return value(); }

    const Error& error() const { return std::get<1>(v_); }

private:
    std::variant<T, Error> v_;
};

}  // namespace geom
