#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
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

int main() {
  // test_find_first_common_node();
  // test_evaluate_rpn();
  // test_tree_is_symmetric();
  // test_flat_preorder_traversal();
  test_find_lowest_common_ancestor();
  return 0;
}
