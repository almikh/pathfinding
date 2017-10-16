#pragma once
#include "path-finder.h"

namespace pf
{
  /// A* path-finder
  template<class T, class Point, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  class PathFinder_AStar: public PathFinder<T, Point> {
  public:
    using HardPtr = std::shared_ptr<PathFinder_AStar>;
    enum Connectivity {
      Four = 4,
      Eight = 8
    };

  protected:
    Connectivity connectivity_ = Eight;
    Matrix<int> label_map_;

  public:
    /// @param connectivity          connectivity of cells (when `Eight` used, diagonal movement is allowed)
    void setConnectivity(Connectivity connectivity) {
      connectivity_ = connectivity;
    }

    /// @param begin                 desired the grid for search
    void setGrid(typename Grid<T>::HardPtr grid) override {
      PathFinder<T, Point>::setGrid(grid);
      label_map_.resize(width_, height_);
    }

    /// @param begin                 starting position 
    /// @param end                   target position
    /// @param include_init_points   need to add start and finish to the path
    bool find(Point start, Point end, bool include_init_points = true) override {
      assert(grid_ != nullptr);
      assert(get_dist_ != nullptr);

      path_.clear();
      total_cost_ = T(0);
      label_map_.clear(0);
      price_.clear(0);

      if (start == end) {
        return true;
      }

      auto comparer = [](const std::pair<Point, T>& left, const std::pair<Point, T>& right) {
        return left.second > right.second;
      };
      std::priority_queue<std::pair<Point, T>, std::vector<std::pair<Point, T>>, decltype(comparer)> queue(comparer);

      int label = 1;
      label_map_(start.x, start.y) = label;
      queue.emplace(Point(start.x, start.y), T(0));

      const float sqrt2 = sqrt(2.0f);

      T temp;
      int nx, ny;
      Point cur, prev(-1, -1);
      for (;;) {
        if (queue.empty()) return false;

        cur = queue.top().first;
        queue.pop();

        if (cur == end) break;

        label = label_map_(cur.x, cur.y);
        for (int i = 0; i < connectivity_; ++i) {
          nx = cur.x + dx[i];
          ny = cur.y + dy[i];
          if (grid_->isTraversableAt(nx, ny) && label_map_(nx, ny) == 0) {
            // get_dist_(nx, ny, cur.x, cur.y) - это фиксированные значения для разных индексов? убрать вычисления?
            temp = price_(cur.x, cur.y) + grid_->getCostAt(nx, ny) + get_dist_(nx, ny, cur.x, cur.y);
            label_map_(nx, ny) = label + 1;

            price_(nx, ny) = temp;
            temp += get_dist_(nx, ny, end.x, end.y);

            queue.emplace(Point(nx, ny), temp);
          }
        }
      }

      // обратная трассировка
      path_.reserve(label + 1);
      if (include_init_points) {
        path_.push_back(end);
      }

      cur = end;
      label = label_map_(cur.x, cur.y) - 1;
      while (cur != start) {
        total_cost_ += grid_->getCostAt(cur.x, cur.y);
        for (int i = 0; i < connectivity_; ++i) {
          nx = cur.x + dx[i];
          ny = cur.y + dy[i];
          if ((start.x == nx && start.y == ny) || (grid_->isTraversableAt(nx, ny) && label_map_(nx, ny) == label)) {
            cur = Point {nx, ny};
            path_.push_back(cur);
            label -= 1;
            break;
          }
        }
      }

      if (!include_init_points) {
        path_.pop_back();
      }

      return true;
    }
  };
}
