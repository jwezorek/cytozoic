#pragma once

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <unordered_set>

namespace cz {

    template <typename T>
    struct vec2 {
        T x;
        T y;

        constexpr bool operator==(const vec2&) const = default;

        constexpr vec2& operator+=(const vec2& rhs) {
            x += rhs.x;
            y += rhs.y;
            return *this;
        }

        constexpr vec2& operator-=(const vec2& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }

        constexpr vec2& operator*=(const T& s) {
            x *= s;
            y *= s;
            return *this;
        }

        constexpr vec2& operator/=(const T& s) {
            x /= s;
            y /= s;
            return *this;
        }
    };

    template <typename T>
    constexpr vec2<T> operator+(vec2<T> lhs, const vec2<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template <typename T>
    constexpr vec2<T> operator-(vec2<T> lhs, const vec2<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template <typename T>
    constexpr vec2<T> operator*(vec2<T> lhs, const T& s) {
        lhs *= s;
        return lhs;
    }

    template <typename T>
    constexpr vec2<T> operator*(const T& s, vec2<T> rhs) {
        rhs *= s;
        return rhs;
    }

    template <typename T>
    constexpr vec2<T> operator/(vec2<T> lhs, const T& s) {
        lhs /= s;
        return lhs;
    }

    template <typename T>
    struct vec2_hash {
        size_t operator()(const vec2<T>& pt) const {
            size_t seed = 0;
            boost::hash_combine(seed, pt.x);
            boost::hash_combine(seed, pt.y);
            return seed;
        }
    };

    template <typename T>
    using vec2_set = std::unordered_set<vec2<T>, vec2_hash<T>>;

    template <typename T, typename U>
    using vec2_map = std::unordered_map<vec2<T>, U, vec2_hash<T>>;

}