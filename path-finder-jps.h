#pragma once
#include "path-finder.h"
#include <map>

namespace pf
{
  /// Path finder using the Jump Point Search algorithm
  template<class T, class Point, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  class PathFinder_JPS: public PathFinder<T, Point> {
  public:
    using HardPtr = std::shared_ptr<PathFinder>;

  protected:
    Matrix<int> label_map_;
    std::map<Point, Point> ancestors_;

    inline int sign(int val) {
      if (val == 0) return 0;
      return (val < 0) ? -1 : 1;
    }

    /// @param current        current node
    /// @param dx             direction for x
    /// @param dy             direction for y
    /// @param begin          starting position
    /// @param end            target position
    virtual Point jump(const Point& current, int dx, int dy, const Point& end) {
      int nx = current.x + dx;
      int ny = current.y + dy;
      for (;;) {
        if (!grid_->isTraversableAt(nx, ny)) return Point(-1, -1);
        if (nx == end.x && ny == end.y) return Point(nx, ny);

        if (dx != 0 && dy != 0) { // diagonal case
          // if ((grid_->isTraversableAt(nx - dx, ny + dy) && !grid_->isTraversableAt(nx - dx, ny)) || 
          //     (grid_->isTraversableAt(nx + dx, ny - dy) && !grid_->isTraversableAt(nx, ny - dy))) {
          //   return Point(nx, ny);
          // }

          // when moving diagonally, must check for vertical/horizontal jump points
          if (jump(Point(nx, ny), dx, 0, end).x >= 0 || jump(Point(nx, ny), 0, dy, end).x >= 0) {
            return Point(nx, ny);
          }
        }
        else if (dx != 0) { // horizontal Forced Neighbor Check
          if ((grid_->isTraversableAt(nx + dx, ny + 1) && !grid_->isTraversableAt(nx, ny + 1)) || 
              (grid_->isTraversableAt(nx + dx, ny - 1) && !grid_->isTraversableAt(nx, ny - 1))) {
            return Point(nx, ny);
          }
        }
        else { // vertical Forced Neighbor Check
          if ((grid_->isTraversableAt(nx + 1, ny + dy) && !grid_->isTraversableAt(nx + 1, ny)) || 
              (grid_->isTraversableAt(nx - 1, ny + dy) && !grid_->isTraversableAt(nx - 1, ny))) {
            return Point(nx, ny);
          }
        }

        nx += dx;
        ny += dy;
      }

      return Point(-1, -1);
    }

    /// @param current        current node
    /// @param begin          starting position
    /// @param end            target position
    virtual std::vector<Point> getSuccesors(const Point& current, const Point& end) {
      int next_x, next_y;
      Point neighbors[8];
      int num_neighbors = 0;
      std::vector<Point> successors;
      if (ancestors_.count(current)) {
        auto parent = ancestors_[current];
        const int x = current.x, y = current.y;
        auto dx = sign(x - parent.x), dy = sign(y - parent.y);
        if (dx != 0 && dy != 0) { // diagonal
          if (grid_->isTraversableAt(x, y + dy)) neighbors[num_neighbors++] = Point(0, dy);
          if (grid_->isTraversableAt(x + dx, y)) neighbors[num_neighbors++] = Point(dx, 0);
          if (grid_->isTraversableAt(x + dx, y + dy)) neighbors[num_neighbors++] = Point(dx, dy);
          // if (grid_->isTraversableAt(x - dx, y + dy)) neighbors[num_neighbors++] = Point(-dx, dy);
          // if (grid_->isTraversableAt(x + dx, y - dy)) neighbors[num_neighbors++] = Point(dx, -dy);
        }
        else if (dx != 0) { // horizontal
          if (grid_->isTraversableAt(x + dx, y)) neighbors[num_neighbors++] = Point(dx, 0);
            if (!grid_->isTraversableAt(x, y + 1)) neighbors[num_neighbors++] = Point(dx, 1);
            if (!grid_->isTraversableAt(x, y - 1)) neighbors[num_neighbors++] = Point(dx, -1);
        }
        else { // vertical
          if (grid_->isTraversableAt(x, y + dy)) neighbors[num_neighbors++] = Point(0, dy);
            if (!grid_->isTraversableAt(x + 1, y)) neighbors[num_neighbors++] = Point(1, dy);
            if (!grid_->isTraversableAt(x - 1, y)) neighbors[num_neighbors++] = Point(-1, dy);
        }
      }
      else {
        for (int i = 0; i < 8; ++i) {
          next_x = current.x + dx[i];
          next_y = current.y + dy[i];
          if (grid_->isTraversableAt(next_x, next_y)) {
            neighbors[num_neighbors++] = Point(dx[i], dy[i]);
          }
        }
      }

      for (int i = 0; i < num_neighbors; ++i) {
        next_x = current.x + neighbors[i].x;
        next_y = current.y + neighbors[i].y;
        auto jump_point = jump(current, neighbors[i].x, neighbors[i].y, end);
        if (jump_point.x >= 0 && !label_map_(jump_point.x, jump_point.y)) { // jump point is found
          successors.push_back(jump_point);
        }
      }

      return successors;
    }

  public:
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
      price_.clear(0);
      label_map_.clear(0);
      total_cost_ = T(0);

      if (start == end) {
        return true;
      }

      auto comparer = [](const std::pair<Point, T>& left, const std::pair<Point, T>& right) {
        return left.second > right.second;
      };
      std::priority_queue<std::pair<Point, T>, std::vector<std::pair<Point, T>>, decltype(comparer)> queue(comparer);

      int label = 1;
      queue.push({Point(start.x, start.y), T(0)});

      T temp;
      int step_x, step_y, step_counter;
      Point cur, temp_pt;
      for (;;) {
        if (queue.empty()) return false;

        cur = queue.top().first;
        queue.pop();

        if (cur == end) break;

        auto succesors = getSuccesors(cur, end);
        for (auto& e : succesors) {
          // compute path cost
          temp = T(0);
          temp_pt = e;
          step_counter = 0;
          step_x = sign(cur.x - e.x);
          step_y = sign(cur.y - e.y);
          while (temp_pt != cur) {
            temp += grid_->getCostAt(temp_pt.x, temp_pt.y);
            temp_pt.x += step_x;
            temp_pt.y += step_y;
            step_counter += 1;
          }

          temp += get_dist_(step_x, step_y, 0, 0) * step_counter;
          temp = price_(cur.x, cur.y) + temp;
          price_(e.x, e.y) = temp;

          ancestors_[e] = cur;
          label_map_(e.x, e.y) = 1;
          queue.push({Point(e.x, e.y), static_cast<T>(temp + get_dist_(e.x, e.y, end.x, end.y))});
        }
      }

      // backtrace
      if (include_init_points) {
        path_.push_back(end);
        total_cost_ = grid_->getCostAt(end.x, end.y);
      }

      cur = end;
      while (cur != start) {
        auto parent = ancestors_[cur];
        auto dx = sign(parent.x - cur.x);
        auto dy = sign(parent.y - cur.y);
        do {
          cur.x += dx;
          cur.y += dy;
          path_.push_back(cur);
          total_cost_ += grid_->getCostAt(cur.x, cur.y);
        } while (cur != parent);
      }
      
      if (!include_init_points) {
        path_.pop_back();
      }

      return true;
    }
  };
}
