#pragma once
#include <limits>
#include <memory>
#include <vector>
#include <cassert>
#include <queue>
#include <cmath>
#include <list>

namespace pf
{
  template<class T>
  struct Grid {
    using HardPtr = std::shared_ptr<Grid<T>>;

    virtual int width() = 0;
    virtual int height() = 0;

    virtual T getCostAt(int x, int y) = 0;
    virtual bool isTraversableAt(int x, int y) = 0;
  };
 
  namespace heuristic
  {
    template<class T>
    T manhattan(int x1, int y1, int x2, int y2) {
      return static_cast<T>(std::abs(x1-x2) + std::abs(y1 - y2));
    }

    template<class T>
    T euclidean(int x1, int y1, int x2, int y2) {
      const int s1 = x1 - x2 , s2 = y1 - y2;
      return std::sqrt(static_cast<T>(s1*s1 + s2*s2));
    }
    
    template<class T>
    T chebyshev(int x1, int y1, int x2, int y2) {
      return static_cast<T>(std::max(std::abs(x1 - x2), std::abs(y1 - y2)));
    }
  }

  template<class T, class Point, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  class PathFinder {
  public:
    using HardPtr = std::shared_ptr<PathFinder>;
    using DistanceFunc = T (*)(int x1, int y1, int x2, int y2);

  // protected
  public: 
    /// Auxiliary class for internal needs
    template<typename S>
    class Matrix { 
      S* data_ = nullptr;
      int width_ = 0, height_ = 0;

      void release() {
        if (data_) {
          delete[] data_;
          data_ = nullptr;
        }

        width_ = height_ = 0;
      }

    public:
      Matrix() = default;
      Matrix(Matrix<S>&& other) = delete;

      Matrix(const Matrix<S>& other):
        data_(new S[other.width_*other.height_]),
        width_(other.width_),
        height_(other.height_)
      {
        memcpy(data_, other.data_, sizeof(S)*width_*height_);
      }

      Matrix(int width, int height, const S& val = 0) {
        resize(width, height);
        clear(val);
      }

      ~Matrix() {
        if (data_) {
          delete[] data_;
        }
      }

      Matrix<S>& operator = (const Matrix<S>& rhs) {
        if (this == &rhs) return *this;

        if (width_*height_ != rhs.width_*rhs.height_) {
          release();
          width_ = rhs.width_;
          height_ = rhs.height_;
          data_ = new S[width_*height_];
        }

        memcpy(data_, rhs.data_, sizeof(S)*width_*height_);
        return *this;
      }

      Matrix<S>& operator = (Matrix<S>&& rhs) = delete;

      void clear(S value) {
        auto ptr = data_;
        const int pixels = width_ * height_;
        for (int i = 0; i < pixels; ++i) {
          *ptr++ = value;
        }
      }

      void resize(int width, int height) {
        if (width_ * height_ != width * height) {
          release();
          data_ = new S[width * height];
        }

        width_ = width;
        height_ = height;
      }

      inline S* data() const {
        return data_;
      }

      inline int height() const {
        return height_;
      }

      inline int width() const {
        return width_;
      }

      inline S& operator () (int x, int y) {
        return data_[x + y*width_];
      }

      inline const S operator () (int x, int y) const {
        return data_[x + y*width_];
      }

      inline S& at(int x, int y) {
        return data_[x + y*width_];
      }

      inline const S at(int x, int y) const {
        return data_[x + y*width_];
      }
    };

  protected:
    const int dx[8] = {-1, 0, 1, 0, -1, 1, 1, -1};
    const int dy[8] = {0, -1, 0, 1, -1, -1, 1, 1};

  protected:
    DistanceFunc get_dist_ = heuristic::euclidean; // used heuristics
    typename Grid<T>::HardPtr grid_; // the grid for search
    uint32_t width_ = 0, height_ = 0; // the grid size
    std::vector<Point> path_; // the found path
    T total_cost_ = 0.0f; // total cost of the found path
    Matrix<T> price_; // cell cost apart from the distance to the target

  public:
    virtual ~PathFinder() = default;

    std::vector<Point> getPath() const {
      return path_;
    }

    T getPathCost() const {
      return total_cost_;
    }

    /// @param begin                 desired heuristic (euclidean by default)
    void setDistanceFunc(typename PathFinder<T, Point>::DistanceFunc func) {
      get_dist_ = func;
    }

    /// @param begin                 desired the grid for search
    virtual void setGrid(typename Grid<T>::HardPtr grid) {
      grid_ = grid;
      width_ = grid_->width();
      height_ = grid_->height();

      price_.resize(width_, height_);
    }

    /// @param begin                 starting position 
    /// @param end                   target position
    /// @param include_init_points   need to add start and finish to the path
    virtual bool find(Point start, Point end, bool include_init_points = true) = 0;
  };
}

#include <path-finder-astar.h>
#include <path-finder-jps.h>
