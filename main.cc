#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <variant>
#include <vector>

template <typename T>
std::ostream &operator<<(std::ostream &strm, const std::vector<T> &v) {
  strm << "[";
  for (const auto &x : v) {
    strm << x << ",";
  }
  strm << "]";
  return strm;
}

template <typename T, size_t n>
std::ostream &operator<<(std::ostream &strm, const std::array<T, n> &v) {
  strm << "[";
  for (const auto &x : v) {
    strm << x << ",";
  }
  strm << "]";
  return strm;
}

struct LinkedListNode {
  std::shared_ptr<LinkedListNode> next;
  int value;
};

std::shared_ptr<LinkedListNode> prepend(int value,
                                        std::shared_ptr<LinkedListNode> tail) {
  return std::make_shared<LinkedListNode>(LinkedListNode{tail, value});
}

int compute_length(std::shared_ptr<LinkedListNode> node) {
  int result = 0;
  while (node != nullptr) {
    result += 1;
    node = node->next;
  }
  return result;
}

std::shared_ptr<LinkedListNode>
find_first_common_node(std::shared_ptr<LinkedListNode> node_a,
                       std::shared_ptr<LinkedListNode> node_b) {
  const int length_a = compute_length(node_a);
  const int length_b = compute_length(node_b);
  for (int i = 0; i < length_b - length_a; ++i) {
    node_b = node_b->next;
  }
  for (int i = 0; i < length_a - length_b; ++i) {
    node_a = node_a->next;
  }
  while (node_a != nullptr) {
    assert(node_b != nullptr);
    if (node_a == node_b) {
      return node_a;
    }
    node_a = node_a->next;
    node_b = node_b->next;
  }
  return nullptr;
}

void test_find_first_common_node() {
  auto common = prepend(1, prepend(0, {}));
  assert(compute_length(common) == 2);
  auto list_a = common;
  auto list_b = common;
  assert(find_first_common_node(list_a, list_b) == common);
  assert(find_first_common_node(list_b, list_a) == common);
  list_a = prepend(2, common);
  assert(compute_length(list_a) == 3);
  assert(find_first_common_node(list_a, list_b) == common);
  assert(find_first_common_node(list_b, list_a) == common);
  list_a = prepend(3, list_a);
  assert(compute_length(list_a) == 4);
  assert(find_first_common_node(list_a, list_b) == common);
  assert(find_first_common_node(list_b, list_a) == common);
  list_b = prepend(4, list_b);
  assert(compute_length(list_b) == 3);
  assert(find_first_common_node(list_a, list_b) == common);
  assert(find_first_common_node(list_b, list_a) == common);
  list_b = prepend(1, {});
  assert(compute_length(list_b) == 1);
  assert(find_first_common_node(list_a, list_b) == nullptr);
  assert(find_first_common_node(list_b, list_a) == nullptr);

  assert(compute_length({}) == 0);
  assert(find_first_common_node(list_a, {}) == nullptr);
  assert(find_first_common_node({}, {}) == nullptr);
}

using RpnToken = std::variant<int, char>;

int evaluate_rpn_op(const int i0, const int i1, const char op) {
  switch (op) {
  case '+':
    return i0 + i1;
  case '-':
    return i0 - i1;
  case '*':
    return i0 * i1;
  case '/':
    return i0 / i1;
  default:
    assert(false);
  }
}

int evaluate_rpn(const std::vector<RpnToken> &tokens) {
  std::stack<int> stack;
  for (const auto &token : tokens) {
    if (const int *i = std::get_if<int>(&token)) {
      stack.push(*i);
    } else if (const char *op = std::get_if<char>(&token)) {
      int i1 = stack.top();
      stack.pop();
      int i0 = stack.top();
      stack.pop();
      stack.push(evaluate_rpn_op(i0, i1, *op));
    }
  }
  assert(stack.size() == 1);
  return stack.top();
}

void test_evaluate_rpn() {
  struct TestCase {
    std::vector<RpnToken> tokens;
    int result;
  };
  std::vector<TestCase> test_cases = {
      {{1, 2, '+'}, 3},
      {{1, 2, 3, '+', '+'}, 6},
      {{1, 2, '+', 3, '-'}, 0},
  };
  for (const auto &[tokens, result_expected] : test_cases) {
    const int result_actual = evaluate_rpn(tokens);
    assert(result_actual == result_expected);
  }
}

struct BinaryTreeNode {
  int value;
  std::shared_ptr<BinaryTreeNode> left;
  std::shared_ptr<BinaryTreeNode> right;
};

std::shared_ptr<BinaryTreeNode>
make_tree(int value, std::shared_ptr<BinaryTreeNode> left,
          std::shared_ptr<BinaryTreeNode> right) {
  return std::make_shared<BinaryTreeNode>(BinaryTreeNode{value, left, right});
}

bool trees_are_mirror_images(std::shared_ptr<BinaryTreeNode> n1,
                             std::shared_ptr<BinaryTreeNode> n2) {
  if (n1 == nullptr && n2 == nullptr) {
    return true;
  }
  if (n1 == nullptr || n2 == nullptr) {
    return false;
  }
  return n1->value == n2->value &&
         trees_are_mirror_images(n1->left, n2->right) &&
         trees_are_mirror_images(n1->right, n2->left);
}

bool tree_is_symmetric(std::shared_ptr<BinaryTreeNode> root) {
  return trees_are_mirror_images(root->left, root->right);
}

void test_tree_is_symmetric() {
  auto *t = &make_tree;
  std::shared_ptr<BinaryTreeNode> tree;
  tree =
      t(0, t(1, t(2, {}, {}), t(3, {}, {})), t(1, t(3, {}, {}), t(2, {}, {})));
  assert(tree_is_symmetric(tree));
  tree =
      t(0, t(1, t(2, {}, {}), t(3, {}, {})), t(2, t(3, {}, {}), t(2, {}, {})));
  assert(!tree_is_symmetric(tree));
  tree = t(0, t(1, t(2, {}, {}), t(3, {}, t(4, {}, {}))),
           t(1, t(3, t(4, {}, {}), {}), t(2, {}, {})));
  assert(tree_is_symmetric(tree));
  tree = t(0, t(1, t(2, {}, {}), t(3, {}, t(4, {}, {}))),
           t(1, t(3, {}, t(4, {}, {})), t(2, {}, {})));
  assert(!tree_is_symmetric(tree));
}

void flat_preorder_traversal(
    std::shared_ptr<BinaryTreeNode> root,
    std::function<bool(std::shared_ptr<BinaryTreeNode>)> action) {
  enum Stage { LEFT, RIGHT, ACTION };
  std::stack<std::pair<std::shared_ptr<BinaryTreeNode>, Stage>> stack;
  auto push_actions = [&stack](auto n) {
    if (n != nullptr) {
      stack.push({n, ACTION});
      stack.push({n, RIGHT});
      stack.push({n, LEFT});
    }
  };
  push_actions(root);
  while (!stack.empty()) {
    const auto [n, stage] = stack.top();
    stack.pop();
    switch (stage) {
    case LEFT:
      push_actions(n->left);
      break;
    case RIGHT:
      push_actions(n->right);
      break;
    case ACTION:
      if (!action(n)) {
        return;
      }
      break;
    }
  }
}

void test_flat_preorder_traversal() {
  auto *t = &make_tree;
  std::ostringstream strm;
  auto track_traversal = [&strm](std::shared_ptr<BinaryTreeNode> n) -> bool {
    strm << n->value << ",";
    return true;
  };

  std::shared_ptr<BinaryTreeNode> tree;
  tree =
      t(0, t(1, t(2, {}, {}), t(3, {}, {})), t(4, t(5, {}, {}), t(6, {}, {})));
  flat_preorder_traversal(tree, track_traversal);
  std::cout << strm.str() << "\n";
}

std::variant<bool, std::shared_ptr<BinaryTreeNode>>
find_lowest_common_ancestor(std::shared_ptr<BinaryTreeNode> n,
                            std::shared_ptr<BinaryTreeNode> n1,
                            std::shared_ptr<BinaryTreeNode> n2) {
  if (n == nullptr) {
    return false;
  }
  auto v_left = find_lowest_common_ancestor(n->left, n1, n2);
  if (auto *lca = std::get_if<1>(&v_left)) {
    return *lca;
  }
  auto v_right = find_lowest_common_ancestor(n->right, n1, n2);
  if (auto *lca = std::get_if<1>(&v_right)) {
    return *lca;
  }
  int n_found = int((n == n1) || (n == n2)) + int(std::get<0>(v_left)) +
                int(std::get<0>(v_right));
  switch (n_found) {
  case 0:
    return false;
  case 1:
    return true;
  case 2:
    return n;
  default:
    assert(false);
  }
}

void test_find_lowest_common_ancestor() {
  auto *t = &make_tree;
  auto root =
      t(0, t(1, t(2, {}, {}), t(3, {}, {})), t(4, t(5, {}, {}), t(6, {}, {})));
  auto res = find_lowest_common_ancestor(root, root->left, root->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root);
  res = find_lowest_common_ancestor(root, root, root->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root);
  res = find_lowest_common_ancestor(root, root, root->right->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root);
  res = find_lowest_common_ancestor(root, root->right, root->right->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root->right);
  res =
      find_lowest_common_ancestor(root, root->right->left, root->right->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root->right);
  res = find_lowest_common_ancestor(root, root->left->left, root->right->right);
  assert(res.index() == 1);
  assert(std::get<1>(res) == root);
}

void sort_up_to_k(std::vector<int> &v, size_t k) {

  std::priority_queue<int, std::vector<int>, std::greater<int>> queue(
      std::greater<int>{});
  int j = 0;
  for (size_t i = 0; i < v.size(); ++i) {
    queue.push(v[i]);
    if (queue.size() == k + 1) {
      v.at(j++) = queue.top();
      queue.pop();
    }
  }
  while (!queue.empty()) {
    v.at(j++) = queue.top();
    queue.pop();
  }
}

template <typename T>
bool vector_equal(const std::vector<T> &v1, const std::vector<T> &v2) {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (size_t i = 0; i < v1.size(); ++i) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

void test_sort_up_to_k() {
  std::vector<int> v{2, 1, 0, 3, 5, 4, 6};
  std::vector<int> w{0, 1, 2, 3, 4, 5, 6};
  sort_up_to_k(v, 3);
  assert(vector_equal(v, w));
}

int find_minimum_cyclically_sorted(const std::vector<int> &v) {
  int i0 = 0;
  int i1 = v.size() / 2;
  int i_min = 0;
  if (v.at(i1) > v.at(i0)) {
    i_min = i0;
    i0 = i1;
    i1 = v.size() - 1;
    // min is in (i1, end) or i0
  } else {
    i_min = i1;
    // min is in (i0, i1]
  }
  while (i0 + 1 < i1) {
    int i2 = (i0 + i1) / 2;
    if (v.at(i2) < v.at(i_min)) {
      i1 = i2;
      i_min = i2;
    } else {
      i0 = i2;
    }
  }
  if (v.at(i1) < v.at(i_min)) {
    i_min = i1;
  }
  return i_min;
}

void test_find_minimum_cyclically_sorted() {
  auto *f = &find_minimum_cyclically_sorted;
  assert(f({1}) == 0);
  assert(f({1, 2}) == 0);
  assert(f({2, 1}) == 1);
  assert(f({1, 2, 3}) == 0);
  assert(f({3, 1, 2}) == 1);
  assert(f({2, 3, 1}) == 2);
  assert(f({1, 2, 3, 4}) == 0);
  assert(f({4, 1, 2, 3}) == 1);
  assert(f({3, 4, 1, 2}) == 2);
  assert(f({2, 3, 4, 1}) == 3);
  assert(f({1, 2, 3, 4, 5}) == 0);
  assert(f({5, 1, 2, 3, 4}) == 1);
  assert(f({4, 5, 1, 2, 3}) == 2);
  assert(f({3, 4, 5, 1, 2}) == 3);
  assert(f({2, 3, 4, 5, 1}) == 4);
}

void reverse(std::vector<int> &v, int begin, int end) {
  int i = begin;
  int j = end - 1;
  while (i < j) {
    std::swap(v[i++], v[j--]);
  }
}

std::vector<int> reversed(const std::vector<int> &v, int begin, int end) {
  auto result = v;
  reverse(result, begin, end);
  return result;
}

void test_reverse() {
  assert(vector_equal(reversed({1, 2, 3, 4}, 0, 4), {4, 3, 2, 1}));
  assert(vector_equal(reversed({1, 2, 3, 4}, 1, 4), {1, 4, 3, 2}));
  assert(vector_equal(reversed({1, 2, 3, 4}, 2, 4), {1, 2, 4, 3}));
  assert(vector_equal(reversed({1, 2, 3, 4}, 0, 3), {3, 2, 1, 4}));
  assert(vector_equal(reversed({1, 2, 3, 4}, 0, 2), {2, 1, 3, 4}));
  assert(vector_equal(reversed({1, 2, 3, 4}, 1, 3), {1, 3, 2, 4}));
}

bool next_permutation(std::vector<int> &v) {
  int i = v.size() - 2;
  while (i >= 0 && v[i] > v[i + 1])
    --i;
  if (i >= 0) {
    int j = i + 1;
    while (j < int(v.size()) && v[j] > v[i]) {
      ++j;
    }
    std::swap(v[i], v[j - 1]);
  }
  reverse(v, i + 1, v.size());
  return i >= 0;
}

void test_next_permutation() {
  auto np = [](const std::vector<int> &v) {
    auto result = v;
    next_permutation(result);
    return result;
  };
  std::vector<std::vector<int>> test_cases = {
      {1, 2, 3, 4, 5}, {1, 2, 3, 5, 4}, {1, 2, 4, 3, 5}, {1, 2, 4, 5, 3},
      {1, 2, 5, 3, 4}, {1, 2, 5, 4, 3}, {1, 3, 2, 4, 5}, {1, 3, 2, 5, 4},
      {1, 3, 4, 2, 5}, {1, 3, 4, 5, 2}, {1, 3, 5, 2, 4}, {1, 3, 5, 4, 2},
  };

  for (size_t i = 0; i + 1 < test_cases.size(); ++i) {
    assert(vector_equal(np(test_cases[i]), test_cases[i + 1]));
  }

  assert(vector_equal(np({5, 4, 3, 2, 1}), {1, 2, 3, 4, 5}));
}

std::shared_ptr<BinaryTreeNode>
find_lowest_common_ancestor_search_tree(std::shared_ptr<BinaryTreeNode> n,
                                        int v1, int v2) {
  if (n == nullptr) {
    return nullptr;
  }
  if (v2 < n->value) {
    return find_lowest_common_ancestor_search_tree(n->left, v1, v2);
  }
  if (v1 > n->value) {
    return find_lowest_common_ancestor_search_tree(n->right, v1, v2);
  }
  return n;
}

void test_find_lowest_common_ancestor_search_tree() {
  auto *t = &make_tree;
  auto root = t(16, t(8, t(4, {}, {}), t(12, {}, {})),
                t(24, t(20, {}, {}), t(28, {}, {})));
  auto res = find_lowest_common_ancestor_search_tree(root, 8, 24);
  assert(res == root);
  res = find_lowest_common_ancestor_search_tree(root, 16, 24);
  assert(res == root);
  res = find_lowest_common_ancestor_search_tree(root, 16, 28);
  assert(res == root);
  res = find_lowest_common_ancestor_search_tree(root, 24, 28);
  assert(res == root->right);
  res = find_lowest_common_ancestor_search_tree(root, 20, 28);
  assert(res == root->right);
  res = find_lowest_common_ancestor_search_tree(root, 4, 28);
  assert(res == root);
}

int count_grid_traversals(int rows, int cols, std::vector<bool> mask) {
  assert(int(mask.size()) == rows * cols);
  std::vector<int> counts(rows * cols, 0);
  counts[0] = 1;
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      const int count = counts[r * cols + c];
      const int i1 = (r + 1) * cols + c;
      const int i2 = r * cols + c + 1;
      if (r + 1 < rows && mask[i1]) {
        counts[i1] += count;
      }
      if (c + 1 < cols && mask[i2]) {
        counts[i2] += count;
      }
    }
  }
  return counts.back();
}

std::optional<int> find_most_common(const std::vector<int> &v, size_t k) {
  std::unordered_map<int, int> counts;
  for (const int x : v) {
    auto pair = counts.find(x);
    if (pair != counts.end()) {
      ++pair->second;
    } else if (counts.size() < k) {
      counts.insert({x, 1});
    } else {
      for (auto pair = counts.begin(); pair != counts.end();) {
        if (pair->second == 1) {
          pair = counts.erase(pair);
        } else {
          --pair->second;
          ++pair;
        }
      }
    }
  }

  std::optional<int> result{};
  for (const auto &pair : counts) {
    if (!result || *result < pair.second) {
      result = pair.first;
    }
  }
  return result;
}

void test_find_most_common() {
  // clang-format off
  std::vector<std::tuple<std::vector<int>, int, int>> test_cases{
      {{1, 1, 1}, 1, 1},
      {{1, 1, 2}, 1, 1},
      {{1, 2, 1}, 1, 1},
      {{2, 1, 1}, 1, 1},
      {{1, 1, 1, 1}, 1, 1},
      {{2, 1, 1, 1}, 1, 1},
      {{1, 2, 1, 1}, 1, 1},
      {{1, 1, 2, 1}, 1, 1},
      {{1, 1, 1, 2}, 1, 1},
      {{1, 1, 1, 1, 1}, 1, 1},
      {{2, 1, 1, 1, 1}, 1, 1},
      {{1, 2, 1, 1, 1}, 1, 1},
      {{1, 1, 2, 1, 1}, 1, 1},
      {{1, 1, 1, 2, 1}, 1, 1},
      {{1, 1, 1, 1, 2}, 1, 1},
      {{2, 2, 1, 1, 1}, 1, 1},
      {{1, 2, 2, 1, 1}, 1, 1},
      {{1, 1, 2, 2, 1}, 1, 1},
      {{1, 1, 1, 2, 2}, 1, 1},
      {{2, 1, 2, 1, 1}, 1, 1},
      {{1, 2, 1, 2, 1}, 1, 1},
      {{1, 1, 2, 1, 2}, 1, 1},
      {{1, 1, 2, 3}, 1, 2},
      {{1, 2, 1, 3}, 1, 2},
      {{2, 1, 1, 3}, 1, 2},
      {{1, 2, 3, 1}, 1, 2},
      {{2, 1, 3, 1}, 1, 2},
      {{2, 3, 1, 1}, 1, 2},
      {{1, 1, 1, 2, 2, 3, 3}, 1, 2},
      {{1, 2, 1, 3, 2, 1, 3}, 1, 2},
  };
  // clang-format on
  for (const auto &[v, expected, k] : test_cases) {
    auto r = find_most_common(v, k);
    assert(r);
    assert(*r == expected);
  }
}

int reverse_digits(int i) {
  const bool negative = i < 0;
  if (negative) {
    i = -i;
  }

  int j = 0;
  while (i > 0) {
    j = 10 * j + (i % 10);
    i = i / 10;
  }

  if (negative) {
    j = -j;
  }
  return j;
}

void test_reverse_digits() {
  assert(reverse_digits(0) == 0);
  assert(reverse_digits(1) == 1);
  assert(reverse_digits(12) == 21);
  assert(reverse_digits(123) == 321);
  assert(reverse_digits(-1) == -1);
  assert(reverse_digits(-12) == -21);
  assert(reverse_digits(-123) == -321);
}

bool sudoku_group_is_duplicate_free(const std::array<char, 81> &grid,
                                    const std::array<int, 9> &group_index) {
  std::array<int, 9> counts{};
  for (const int i : group_index) {
    if (grid.at(i) > 0 && (++counts.at(grid.at(i) - 1)) > 1) {
      return false;
    }
  }
  return true;
}

constexpr auto make_sudoku_groups() {
  std::array<std::array<int, 9>, 27> groups{};
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      groups[i][j] = 9 * i + j;
      groups[i + 9][j] = 9 * j + i;
      groups[i + 18][j] = 9 * ((i / 3) * 3 + (j / 3)) + ((i % 3) * 3 + (j % 3));
    }
  }
  return groups;
}

bool sudoku_is_valid(const std::array<char, 81> &grid) {
  constexpr auto groups = make_sudoku_groups();
  for (const auto &group : groups) {
    if (!sudoku_group_is_duplicate_free(grid, group)) {
      return false;
    }
  }
  return true;
}

void test_sudoku_is_valid() {
  const auto groups = make_sudoku_groups();
  for (const auto &group : groups) {
    std::cout << std::vector(group.begin(), group.end()) << "\n";
  }
  std::array<char, 81> grid{};
  assert(sudoku_is_valid(grid));
  grid[9 * 0 + 0] = 1;
  assert(sudoku_is_valid(grid));
  grid[9 * 0 + 1] = 1;
  assert(!sudoku_is_valid(grid));
  grid[9 * 0 + 1] = 2;
  assert(sudoku_is_valid(grid));
  grid[9 * 0 + 3] = 1;
  assert(!sudoku_is_valid(grid));
  grid[9 * 0 + 3] = 3;
  assert(sudoku_is_valid(grid));
  grid[9 * 1 + 0] = 1;
  assert(!sudoku_is_valid(grid));
  grid[9 * 1 + 0] = 4;
  assert(sudoku_is_valid(grid));
}

void flood_fill(int rows, int cols, std::vector<bool> &image, int r, int c) {
  assert(int(image.size()) == rows * cols);
  constexpr std::array<std::array<int, 2>, 4> deltas{
      {{1, 0}, {0, 1}, {-1, 0}, {0, -1}}};
  const bool color = image[r * rows + c];
  std::queue<std::array<int, 2>> queue;
  image[r * cols + c] = !color;
  queue.push({r, c});
  while (!queue.empty()) {
    const auto [r, c] = queue.front();
    queue.pop();
    for (const auto &[dr, dc] : deltas) {
      const int r1 = r + dr;
      const int c1 = c + dc;
      if (r1 >= 0 && r1 < rows && c1 >= 0 && c1 < cols &&
          image[r1 * cols + c1] == color) {
        image[r1 * cols + c1] = !color;
        queue.push({r1, c1});
      }
    }
  }
}

void test_flood_fill() {
  // clang-format off
  std::vector<bool> image = {
      0,0,1,0,0,0,
      0,0,1,0,0,0,
      0,0,1,0,0,0,
      1,1,1,0,0,0,
      0,0,0,0,0,0,
  };
  std::vector<bool> image_expected = {
      1,1,1,0,0,0,
      1,1,1,0,0,0,
      1,1,1,0,0,0,
      1,1,1,0,0,0,
      0,0,0,0,0,0,
  };
  // clang-format on
  flood_fill(5, 6, image, 0, 1);
  assert(vector_equal(image, image_expected));
}

void reverse_word_order(std::string &str) {
  std::reverse(str.begin(), str.end());
  auto it1 = str.begin();
  while (it1 != str.end()) {
    while (it1 != str.end() && *it1 == ' ')
      ++it1;
    auto it2 = std::next(it1);
    while (it2 != str.end() && *it2 != ' ')
      ++it2;
    std::reverse(it1, it2);
    it1 = it2;
  }
}

void test_reverse_word_order() {
  auto rev = [](const std::string &str) {
    std::string res = str;
    reverse_word_order(res);
    return res;
  };
  assert(rev("Alice likes Bob") == "Bob likes Alice");
}

std::shared_ptr<LinkedListNode>
remove_kth_last(const std::shared_ptr<LinkedListNode> n, int k) {
  auto q = n;
  auto p = n;
  for (int i = 0; i < k + 1; ++i) {
    if (q == nullptr) {
      return nullptr;
    }
    q = q->next;
  }
  if (q == nullptr) {
    return p->next;
  }
  while (q->next != nullptr) {
    q = q->next;
    p = p->next;
  }
  p->next = p->next->next;
  return n;
}

std::vector<int> to_vector(std::shared_ptr<LinkedListNode> n) {
  std::vector<int> result;
  while (n != nullptr) {
    result.push_back(n->value);
    n = n->next;
  }
  return result;
}

void test_remove_kth_last() {
  auto p = &prepend;
  assert(remove_kth_last(nullptr, 0) == nullptr);
  assert(remove_kth_last(nullptr, 1) == nullptr);
  assert(remove_kth_last(p(1, nullptr), 0) == nullptr);
  assert(remove_kth_last(p(1, nullptr), 1) == nullptr);
  assert(remove_kth_last(p(1, nullptr), 2) == nullptr);
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, nullptr)), 0)), {1}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, nullptr)), 1)), {2}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, nullptr)), 2)), {}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, p(3, nullptr))), 0)),
                      {1, 2}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, p(3, nullptr))), 1)),
                      {1, 3}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, p(3, nullptr))), 2)),
                      {2, 3}));
  assert(vector_equal(to_vector(remove_kth_last(p(1, p(2, p(3, nullptr))), 3)),
                      {}));
}

template <typename T> class CircularQueue {
public:
  CircularQueue(const size_t capacity)
      : data_(capacity, {}), begin_{}, end_{}, empty_{true} {
    if (capacity == 0) {
      throw std::logic_error("Zero capacity");
    }
  }

  size_t size() const {
    if (end_ > begin_)
      return end_ - begin_;
    if (end_ < begin_)
      return end_ + data_.size() - begin_;
    if (empty_)
      return 0;
    return data_.size();
  }

  size_t capacity() const { return data_.size(); }

  void push(T &&x) {
    if (!empty_ && end_ == begin_) {
      std::vector<T> new_data(data_.capacity() * 2);
      auto it =
          std::copy(data_.begin() + begin_, data_.end(), new_data.begin());
      std::copy(data_.begin(), data_.begin() + begin_, it);
      begin_ = 0;
      end_ = data_.size();
      data_ = std::move(new_data);
    }
    empty_ = false;
    data_[end_++] = std::move(x);
    if (end_ == data_.size()) {
      end_ = 0;
    }
  }

  T pop() {
    if (empty_) {
      throw std::logic_error("Empty");
    }
    auto result = std::move(data_[begin_++]);
    if (begin_ == data_.size()) {
      begin_ = 0;
    }
    empty_ = (begin_ == end_);
    return result;
  }

private:
  std::vector<T> data_;
  size_t begin_;
  size_t end_;
  bool empty_;
};

template <typename Func> void assert_throws(Func func) {
  bool caught = false;
  try {
    func();
  } catch (std::logic_error &e) {
    caught = true;
  }

  assert(caught);
}

void test_circular_queue() {
  {
    CircularQueue<int> queue(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }

  {
    CircularQueue<int> queue(1);
    queue.push(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    assert(queue.pop() == 1);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }

  {
    CircularQueue<int> queue(1);
    queue.push(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    queue.push(2);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    assert(queue.pop() == 1);
    assert(queue.size() == 1);
    assert(queue.pop() == 2);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }

  {
    CircularQueue<int> queue(1);
    queue.push(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    queue.push(2);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    queue.push(3);
    assert(queue.capacity() == 4);
    assert(queue.size() == 3);
    assert(queue.pop() == 1);
    assert(queue.size() == 2);
    assert(queue.pop() == 2);
    assert(queue.size() == 1);
    assert(queue.pop() == 3);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }

  {
    CircularQueue<int> queue(1);
    queue.push(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    assert(queue.pop() == 1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 0);
    queue.push(2);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    queue.push(3);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    assert(queue.pop() == 2);
    assert(queue.size() == 1);
    assert(queue.pop() == 3);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }

  {
    CircularQueue<int> queue(1);
    queue.push(1);
    assert(queue.capacity() == 1);
    assert(queue.size() == 1);
    queue.push(2);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    assert(queue.pop() == 1);
    assert(queue.capacity() == 2);
    assert(queue.size() == 1);
    assert(queue.pop() == 2);
    assert(queue.capacity() == 2);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
    queue.push(3);
    assert(queue.capacity() == 2);
    assert(queue.size() == 1);
    queue.push(4);
    assert(queue.capacity() == 2);
    assert(queue.size() == 2);
    assert(queue.pop() == 3);
    assert(queue.size() == 1);
    assert(queue.pop() == 4);
    assert(queue.size() == 0);
    assert_throws([&queue]() { queue.pop(); });
  }
}

std::shared_ptr<BinaryTreeNode>
reconstruct_tree(const int *const inorder_begin, const int *const inorder_end,
                 const int *const preorder_begin,
                 const int *const preorder_end) {
  if ((inorder_end < inorder_begin) ||
      ((inorder_end - inorder_begin) != (preorder_end - preorder_begin))) {
    throw std::logic_error("Invalid ranges");
  }
  if (inorder_begin == inorder_end) {
    return nullptr;
  }

  const int root_value = *preorder_begin;
  auto it = std::find(inorder_begin, inorder_end, root_value);
  if (it == inorder_end) {
    throw std::logic_error("Incompatible traversals");
  }
  const size_t n_left = it - inorder_begin;
  return make_tree(root_value,
                   reconstruct_tree(inorder_begin, inorder_begin + n_left,
                                    preorder_begin + 1,
                                    preorder_begin + 1 + n_left),
                   reconstruct_tree(inorder_begin + n_left + 1, inorder_end,
                                    preorder_begin + n_left + 1, preorder_end));
}

std::optional<int16_t> find_one_missing(const std::vector<uint16_t> &v) {
  std::array<int, 256> counts{};
  const uint16_t mask = 0xFF;
  for (const uint16_t x : v) {
    ++counts.at(x & mask);
  }
  size_t i = 0;
  for (i = 0; i < counts.size(); ++i) {
    if (counts[i] <= 0xFF)
      break;
  }
  if (i == counts.size()) {
    return {};
  }
  counts.fill(0);
  for (const uint16_t x : v) {
    if ((x & mask) == i) {
      ++counts.at((x >> 8) & mask);
    }
  }
  size_t j = 0;
  for (j = 0; j < counts.size(); ++j) {
    if (counts[j] == 0)
      break;
  }
  assert(j < counts.size());
  return (j << 8) | i;
}

void test_find_one_missing() {
  std::vector<uint16_t> v(0xFFFF, 0);
  std::iota(v.begin(), v.end(), 0);
  std::mt19937 rng;
  std::shuffle(v.begin(), v.end(), rng);
  const size_t i = std::uniform_int_distribution<size_t>(0, v.size() - 1)(rng);
  const uint16_t missing = v[i];
  v.erase(v.begin() + i);
  assert(find_one_missing(v).value() == missing);
}

std::array<size_t, 2> find_closest_repetition(const std::vector<int> &v) {
  std::unordered_map<size_t, int> last_seen;
  size_t min_distance = -1;
  std::array<size_t, 2> result{};
  for (size_t i = 0; i < v.size(); ++i) {
    const auto &x = v[i];
    auto it = last_seen.insert({x, i}).first;
    const size_t distance = i - it->second;
    if (distance > 0 && distance < min_distance) {
      min_distance = distance;
      result[0] = it->second;
      result[1] = i;
    }
    it->second = i;
  }
  return result;
}

void test_find_closest_repetition() {
  auto test_case = [](const std::vector<int> &v,
                      std::array<size_t, 2> expected) {
    auto res = find_closest_repetition(v);
    assert(res[0] == expected[0]);
    assert(res[1] == expected[1]);
  };
  test_case({1, 2, 3, 1, 4, 5}, {0, 3});
  test_case({1, 2, 2, 1, 4, 5}, {1, 2});
}

int maximum_concurrent_events(const std::vector<std::array<int, 2>> &events) {
  std::vector<std::pair<int, bool>> end_points;
  end_points.reserve(events.size() * 2);
  for (const auto &[begin, end] : events) {
    end_points.push_back({begin, false});
    end_points.push_back({end, true});
  }
  std::sort(end_points.begin(), end_points.end());
  int n_concurrent = 0;
  int max_n_concurrent = 0;
  for (const auto &[time, is_end] : end_points) {
    if (is_end) {
      --n_concurrent;
    } else {
      ++n_concurrent;
    }
    if (n_concurrent > max_n_concurrent) {
      max_n_concurrent = n_concurrent;
    }
  }
  return max_n_concurrent;
}

void test_maximum_concurrent_events() {
  assert(maximum_concurrent_events({}) == 0);
  assert(maximum_concurrent_events({{0, 1}}) == 1);
  assert(maximum_concurrent_events({{0, 1}, {2, 3}}) == 1);
  assert(maximum_concurrent_events({{0, 1}, {1, 2}, {2, 3}}) == 2);
  assert(maximum_concurrent_events({{0, 2}, {1, 2}, {3, 4}}) == 2);
  assert(maximum_concurrent_events({{0, 1}, {2, 3}, {3, 4}}) == 2);
  assert(maximum_concurrent_events({{0, 1}, {2, 3}, {4, 5}}) == 1);
  assert(maximum_concurrent_events({{0, 2}, {1, 2}, {2, 3}, {3, 4}}) == 3);
}

class MostCommonTokens {
public:
  MostCommonTokens(size_t max_size) : max_size_(max_size) {}

  void insert(const std::string &token) {
    auto token_and_count = token_to_count_.insert({token, 0}).first;
    ++token_and_count->second;
    auto token_and_it = token_to_iterator_.find(token);
    if (token_and_it != token_to_iterator_.end()) {
      count_to_token_.erase(token_and_it->second);
      token_to_iterator_.erase(token_and_it);
    }
    token_to_iterator_.insert(
        {token, count_to_token_.insert(
                    {token_and_count->second, token_and_count->first})});
    if (count_to_token_.size() > max_size_) {
      token_to_iterator_.erase(count_to_token_.begin()->second);
      count_to_token_.erase(count_to_token_.begin());
    }
  }

  const std::multimap<int, std::string> &count_to_token() {
    return count_to_token_;
  }

private:
  size_t max_size_;
  std::unordered_map<std::string, int> token_to_count_;
  std::multimap<int, std::string> count_to_token_;
  std::unordered_map<std::string, std::map<int, std::string>::iterator>
      token_to_iterator_;
};

struct KnapsackSolution {
  int value;
  std::vector<size_t> objects;
};
auto knapsack(const std::vector<int> &w, const std::vector<int> &v,
              const int max_w) {
  std::map<int, KnapsackSolution> wtos{{0, {0, {}}}};
  for (size_t i = 0; i < w.size(); ++i) {
    assert(w[i] > 0);
    for (auto it = wtos.begin(); it != wtos.end(); ++it) {
      auto &[w0, s0] = *it;
      if (!s0.objects.empty() && s0.objects.back() == i) {
        continue;
      }
      const int w1 = w0 + w[i];
      if (w1 > max_w) {
        continue;
      }
      const int v1 = s0.value + v[i];
      auto ub = wtos.upper_bound(w1);
      assert(ub != wtos.begin());
      auto it1 = std::prev(ub);
      assert(it1->first <= w1);
      if (it1->second.value > v1) {
        continue;
      }
      if (it1->first != w1) {
        const auto [it2, inserted] = wtos.insert({w1, s0});
        it1 = it2;
        assert(inserted);
      }
      it1->second.value = v1;
      it1->second.objects.push_back(i);
      ++it1;
      while (it1 != wtos.end() && it1->second.value < v1)
        it1 = wtos.erase(it1);
    }
  }
  return wtos;
}

void test_knapsack() {
  const std::vector<std::vector<int>> ww = {
      {8, 4, 5},
      {8, 4, 2, 3},
      {8, 4, 2, 3, 3},
  };

  const std::vector<std::vector<int>> vv = {
      {32, 15, 19},
      {32, 15, 7, 11},
      {32, 15, 7, 11, 10},
  };

  for (size_t i = 0; i < ww.size(); ++i) {
    auto wtos = knapsack(ww[i], vv[i], 1024);
    std::cout << "w: " << ww[i] << "\nv: " << vv[i] << "\nSolution:\n";
    for (const auto &[w, s] : wtos) {
      std::cout << "w: " << w << " v: " << s.value << " obj: " << s.objects
                << "\n";
    }
    std::cout << "\n";
  }
}

std::vector<std::array<size_t, 2>> max_water(const std::vector<int> &a) {
  std::vector<size_t> hi = {0};
  for (size_t i = 1; i < a.size(); ++i) {
    if (a[i] > a[hi.back()]) {
      hi.push_back(i);
    }
  }

  int max_area = 0;
  std::vector<std::array<size_t, 2>> i_max_area;

  int max_a = 0;
  for (size_t i = a.size() - 1; i < a.size(); --i) {
    if (a[i] <= max_a)
      continue;
    for (const auto &j : hi) {
      if (j >= i) {
        break;
      }
      const int area = (i - j) * std::min(a[i], a[j]);
      if (area > max_area) {
        i_max_area.clear();
      }
      if (area >= max_area) {
        max_area = area;
        i_max_area.push_back({j, i});
      }
    }
  }
  return i_max_area;
}

std::vector<std::array<size_t, 2>>
max_water_slow_and_correct(const std::vector<int> &a) {
  int max_area = 0;
  std::vector<std::array<size_t, 2>> i_max_area;
  for (size_t i = 0; i < a.size(); ++i) {
    for (size_t j = i + 1; j < a.size(); ++j) {
      const int area = (j - i) * std::min(a[i], a[j]);
      if (area > max_area) {
        i_max_area.clear();
      }
      if (area >= max_area) {
        max_area = area;
        i_max_area.push_back({i, j});
      }
    }
  }
  return i_max_area;
}

void test_max_water() {
  std::mt19937 rng;
  for (size_t it = 0; it < 128; ++it) {
    std::uniform_int_distribution uniform_int(0, 16);
    std::vector<int> a(32);
    for (size_t i = 0; i < a.size(); ++i) {
      a[i] = uniform_int(rng);
    }
    auto v1 = max_water(a);
    auto v2 = max_water_slow_and_correct(a);
    std::cout << "a:  " << a << "\n";
    std::cout << "v1: " << v1 << "\n";
    std::cout << "v2: " << v2 << "\n";
  }
}

auto largest_building(const std::vector<int> &h) {
  std::stack<std::pair<size_t, int>> open;
  int max_area = 0;
  std::vector<std::array<size_t, 3>> i_max_area;
  for (size_t i = 0; i < h.size() + 1; ++i) {
    const int h_i = i < h.size() ? h[i] : 0;
    size_t i_start = i;
    while (!open.empty() && open.top().second >= h_i) {
      i_start = open.top().first;
      const int area = (i - i_start) * open.top().second;
      if (area > max_area) {
        i_max_area.clear();
      }
      if (area >= max_area) {
        max_area = area;
        i_max_area.push_back({i_start, i, size_t(open.top().second)});
      }
      open.pop();
    }
    open.push({i_start, h_i});
  }
  return i_max_area;
}

auto largest_building_slow_and_correct(const std::vector<int> &h) {
  int max_area = 0;
  std::vector<std::array<size_t, 3>> i_max_area;
  for (size_t i = 0; i < h.size(); ++i) {
    for (size_t j = i; j < h.size(); ++j) {
      int min_h = h[i];
      for (size_t k = i; k <= j; ++k) {
        min_h = std::min(min_h, h[k]);
      }
      const int area = (j - i + 1) * min_h;
      if (area > max_area) {
        i_max_area.clear();
      }
      if (area >= max_area) {
        max_area = area;
        i_max_area.push_back({i, j + 1, size_t(min_h)});
      }
    }
  }
  return i_max_area;
}

void test_largest_building() {
  std::mt19937 rng;
  for (size_t it = 0; it < 128; ++it) {
    std::uniform_int_distribution uniform_int(0, 16);
    std::vector<int> h(32);
    for (size_t i = 0; i < h.size(); ++i) {
      h[i] = uniform_int(rng);
    }
    auto i1 = largest_building(h);
    auto i2 = largest_building_slow_and_correct(h);
    std::cout << "h:  " << h << "\n";
    std::cout << "i1: " << i1 << "\n";
    std::cout << "i2: " << i2 << "\n";
  }
}

void mark_interior(int rows, int cols, std::vector<char> &img) {
  using Ix = std::array<int, 2>;
  const std::vector<Ix> displ = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  auto at = [&img, cols](const Ix &ix) -> char & {
    return img.at(ix[0] * cols + ix[1]);
  };
  std::queue<Ix> queue;
  auto maybe_push = [&at, &img, &queue, cols](const Ix &ix) -> void {
    char &pixel = at(ix);
    if (pixel == 'W') {
      pixel = 'V';
      queue.push(ix);
    }
  };
  for (int r = 0; r < rows; ++r) {
    maybe_push({r, 0});
    maybe_push({r, cols - 1});
  }
  for (int c = 0; c < cols; ++c) {
    maybe_push({0, c});
    maybe_push({rows - 1, c});
  }

  while (!queue.empty()) {
    const Ix ix = queue.front();
    queue.pop();
    for (const auto &[dr, dc] : displ) {
      const Ix ix1{ix[0] + dr, ix[1] + dc};
      if (ix1[0] >= 0 && ix1[0] < rows && ix1[1] >= 0 && ix1[1] < cols) {
        maybe_push(ix1);
      }
    }
  }

  for (char &pixel : img) {
    if (pixel == 'V') {
      pixel = 'W';
    } else if (pixel == 'W') {
      pixel = 'B';
    }
  }
}

void test_mark_interior() {
  // clang-format off
  std::vector<char> img = {
    'W','W','B','W','W','W',
    'W','B','B','W','B','W',
    'B','W','W','B','W','B',
    'W','B','B','W','B','W',
    'W','W','B','W','W','W',
  };
  // clang-format on
  const int rows = 5;
  const int cols = 6;
  mark_interior(rows, cols, img);
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      std::cout << img.at(r * cols + c) << ",";
    }
    std::cout << "\n";
  }
}

bool sequence_appears_in_matrix_recursive(const std::vector<int> &seq,
                                          const size_t seq_start,
                                          const size_t row, const size_t col,
                                          const size_t rows, const size_t cols,
                                          const std::vector<int> &matrix) {
  const std::vector<std::array<int, 2>> displ = {
      {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  if (seq_start == seq.size()) {
    return true;
  }
  if (matrix.at(row * cols + col) == seq.at(seq_start)) {
    for (const auto [dr, dc] : displ) {
      const size_t row1 = row + dr;
      const size_t col1 = col + dc;
      if (row1 < rows && col1 < cols &&
          sequence_appears_in_matrix_recursive(seq, seq_start + 1, row1, col1,
                                               rows, cols, matrix)) {
        return true;
      }
    }
  }
  return false;
}

bool sequence_appears_in_matrix(const std::vector<int> &seq, const size_t rows,
                                const size_t cols,
                                const std::vector<int> &matrix) {
  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      if (sequence_appears_in_matrix_recursive(seq, 0, row, col, rows, cols,
                                               matrix)) {
        return true;
      }
    }
  }
  return false;
}

void test_sequence_appears_in_matrix() {
  // clang-format off
  const std::vector<int> matrix = {
    0, 1, 2, 3, 4, 5,
    6, 7, 8, 9, 0, 1,
    2, 3, 4, 5, 6, 7,
    8, 9, 0, 1, 2, 3,
  };
  // clang-format on
  const int rows = 4;
  const int cols = 6;

  assert(sequence_appears_in_matrix({0, 1, 7, 8, 4}, rows, cols, matrix));
  assert(sequence_appears_in_matrix({8, 9, 3, 4, 0}, rows, cols, matrix));
  assert(!sequence_appears_in_matrix({8, 9, 3, 4, 9}, rows, cols, matrix));
}

int main() {
  // test_find_first_common_node();
  // test_evaluate_rpn();
  // test_tree_is_symmetric();
  // test_flat_preorder_traversal();
  // test_sort_up_to_k();
  // test_find_minimum_cyclically_sorted();
  // test_find_lowest_common_ancestor();
  // test_flat_tree_traverse();
  // test_reverse();
  // test_next_permutation();
  // test_find_minimum_cyclically_sorted();
  // test_find_lowest_common_ancestor_search_tree();
  // test_find_most_common();
  // test_reverse_digits();
  // test_sudoku_is_valid();
  // test_flood_fill();
  // test_reverse_word_order();
  // test_remove_kth_last();
  // test_circular_queue();
  // test_find_one_missing();
  // test_find_closest_repetition();
  // test_maximum_concurrent_events();
  // test_knapsack();
  // test_max_water();
  // test_largest_building();
  // test_mark_interior();
  test_sequence_appears_in_matrix();

  return 0;
}
