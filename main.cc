#include <cassert>
#include <memory>

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

int main() {
  test_find_first_common_node();
  return 0;
}
