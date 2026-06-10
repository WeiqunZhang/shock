#ifndef FAB_H_
#define FAB_H_

#include "Extension.H"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <utility>

#if defined(DEBUG)
#include <iostream>
#endif

template <typename T>
class FAB
{
public:
    using value_type = T;

    FAB () = default;

    FAB (std::array<int,3> const& lo, // inclusive lower bounds
         std::array<int,3> const& hi, // inclusive upper bounds
         int n = 1);                  // number of components

    ~FAB ();

    FAB (FAB<T> const&) = delete;
    FAB<T>& operator= (FAB<T> const&) = delete;

    FAB (FAB<T>&& rhs) noexcept;
    FAB<T>& operator= (FAB<T>&& rhs) noexcept;

    void clear () noexcept;

    [[nodiscard]] int nComp () const noexcept { return ncomp; }

    [[nodiscard]] std::array<int,3> const& lo () const noexcept { return begin; }

    [[nodiscard]] std::array<int,3> hi () const noexcept
    {
        return { end[0] - 1, end[1] - 1, end[2] - 1 };
    }

    [[nodiscard]] T* data () noexcept { return p; }
    [[nodiscard]] T const* data () const noexcept { return p; }

    [[nodiscard]] FORCE_INLINE
    T& operator() (int i, int j, int k, int n = 0) noexcept
    {
#if defined(DEBUG)
        bound_check(i, j, k, n);
#endif
        return p[index(i, j, k, n)];
    }

    [[nodiscard]] FORCE_INLINE
    T const& operator() (int i, int j, int k, int n = 0) const noexcept
    {
#if defined(DEBUG)
        bound_check(i, j, k, n);
#endif
        return p[index(i, j, k, n)];
    }

private:
#if defined(DEBUG)
    FORCE_INLINE
    void bound_check (int i, int j, int k, int n) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1] ||
            k < begin[2] || k >= end[2] ||
            n < 0          || n >= ncomp)
        {
            std::cerr << "FAB index (" << i << "," << j << "," << k << "," << n
                      << ") out of bounds ("
                      << begin[0] << ":" << end[0] - 1 << ","
                      << begin[1] << ":" << end[1] - 1 << ","
                      << begin[2] << ":" << end[2] - 1 << ","
                      << "0:" << ncomp - 1 << ")\n";
            std::abort();
        }
    }
#endif

    [[nodiscard]] FORCE_INLINE
    std::int64_t index (int i, int j, int k, int n) const noexcept
    {
        return static_cast<std::int64_t>((i - begin[0])
                                       + (j - begin[1]) * stride[0]
                                       + (k - begin[2]) * stride[1]
                                       +  n             * stride[2]);
    }

    static void set_stride (std::array<std::int64_t,3>& stride,
                            std::array<int,3> const& begin,
                            std::array<int,3> const& end) noexcept
    {
        std::int64_t s = 1;
        for (int d = 0; d < 3; ++d) {
            s *= static_cast<std::int64_t>(end[d] - begin[d]);
            stride[d] = s;
        }
    }

    T* RESTRICT p = nullptr;
    std::array<std::int64_t,3> stride{};
    std::array<int,3> begin{}; // inclusive lower bounds
    std::array<int,3> end{};   // exclusive upper bounds
    int ncomp = 0;
    bool owner = false;
};

template <typename T>
FAB<T>::FAB (std::array<int,3> const& lo,
             std::array<int,3> const& hi,
             int n)
    : begin(lo)
    , end{ hi[0] + 1, hi[1] + 1, hi[2] + 1 }
    , ncomp(n)
    , owner(true)
{
    set_stride(stride, begin, end);

    if (ncomp > 0 && stride[0] > 0 && stride[1] > 0 && stride[2] > 0) {
        std::int64_t const nelts = stride[2] * static_cast<std::int64_t>(ncomp);
        p = static_cast<T*>(std::malloc(static_cast<std::size_t>(nelts) * sizeof(T)));
#if defined(DEBUG)
        if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
            T const qnan = std::numeric_limits<T>::quiet_NaN();
            for (std::int64_t i = 0; i < nelts; ++i) {
                p[i] = qnan;
            }
        }
#endif
    }
}

template <typename T>
void FAB<T>::clear () noexcept
{
    if (owner && p != nullptr) {
        std::free(p);
    }
    p = nullptr;
    stride = {};
    begin = {};
    end = {};
    ncomp = 0;
    owner = false;
}

template <typename T>
FAB<T>::~FAB ()
{
    clear();
}

template <typename T>
FAB<T>::FAB (FAB<T>&& rhs) noexcept
    : p(rhs.p)
    , stride(rhs.stride)
    , begin(rhs.begin)
    , end(rhs.end)
    , ncomp(rhs.ncomp)
    , owner(rhs.owner)
{
    rhs.p = nullptr;
    rhs.owner = false;
}

template <typename T>
FAB<T>& FAB<T>::operator= (FAB<T>&& rhs) noexcept
{
    if (this != &rhs) {
        clear();

        p = rhs.p;
        stride = rhs.stride;
        begin = rhs.begin;
        end = rhs.end;
        ncomp = rhs.ncomp;
        owner = rhs.owner;

        rhs.p = nullptr;
        rhs.owner = false;
    }
    return *this;
}

#endif
