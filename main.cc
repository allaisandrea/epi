#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <stack>
#include <variant>
#include <vector>

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

bool vector_equal(const std::vector<int> &v1, const std::vector<int> &v2) {
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
  test_next_permutation();
  return 0;
}
