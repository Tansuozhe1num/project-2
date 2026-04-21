/** * collision_world.c -- OpenCilk Parallelized Quadtree Collision Detection
 * Project 2 - Part 1 & Part 2
 **/

#include "./collision_world.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

// 引入 OpenCilk 核心头文件
#include <cilk/cilk.h>

#include "./intersection_detection.h"
#include "./intersection_event_list.h"
#include "./line.h"

// =============================================================================
// 四叉树并行化支持模块
// =============================================================================

#define QT_MAX_DEPTH 8
#define QT_MAX_LINES 16

typedef struct QuadtreeNode {
  double x_min, y_min, x_max, y_max;
  Line** lines;
  int line_count;
  int capacity;
  struct QuadtreeNode* children[4];
  bool is_leaf;
} QuadtreeNode;

QuadtreeNode* QuadtreeNode_new(double x_min, double y_min, double x_max, double y_max) {
  QuadtreeNode* node = malloc(sizeof(QuadtreeNode));
  node->x_min = x_min; node->y_min = y_min;
  node->x_max = x_max; node->y_max = y_max;
  node->capacity = QT_MAX_LINES;
  node->lines = malloc(sizeof(Line*) * node->capacity);
  node->line_count = 0;
  node->is_leaf = true;
  for (int i = 0; i < 4; i++) node->children[i] = NULL;
  return node;
}

void QuadtreeNode_delete(QuadtreeNode* node) {
  if (!node->is_leaf) {
    for (int i = 0; i < 4; i++) {
      if (node->children[i]) QuadtreeNode_delete(node->children[i]);
    }
  }
  free(node->lines);
  free(node);
}

bool lineOverlapsBox(Line* l, double dt, double xmin, double ymin, double xmax, double ymax) {
  double l_xmin = fmin(l->p1.x, l->p2.x);
  double l_xmax = fmax(l->p1.x, l->p2.x);
  double l_ymin = fmin(l->p1.y, l->p2.y);
  double l_ymax = fmax(l->p1.y, l->p2.y);

  double nx1 = l->p1.x + l->velocity.x * dt;
  double nx2 = l->p2.x + l->velocity.x * dt;
  double ny1 = l->p1.y + l->velocity.y * dt;
  double ny2 = l->p2.y + l->velocity.y * dt;

  l_xmin = fmin(l_xmin, fmin(nx1, nx2));
  l_xmax = fmax(l_xmax, fmax(nx1, nx2));
  l_ymin = fmin(l_ymin, fmin(ny1, ny2));
  l_ymax = fmax(l_ymax, fmax(ny1, ny2));

  return !(l_xmax < xmin || l_xmin > xmax || l_ymax < ymin || l_ymin > ymax);
}

void Quadtree_insert(QuadtreeNode* node, Line* l, int depth, double dt) {
  if (!lineOverlapsBox(l, dt, node->x_min, node->y_min, node->x_max, node->y_max)) return;

  if (node->is_leaf) {
    if (node->line_count < node->capacity || depth >= QT_MAX_DEPTH) {
      if (node->line_count == node->capacity) {
        node->capacity *= 2;
        node->lines = realloc(node->lines, sizeof(Line*) * node->capacity);
      }
      node->lines[node->line_count++] = l;
      return;
    }
    node->is_leaf = false;
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;

    node->children[0] = QuadtreeNode_new(node->x_min, node->y_min, mid_x, mid_y);
    node->children[1] = QuadtreeNode_new(mid_x, node->y_min, node->x_max, mid_y);
    node->children[2] = QuadtreeNode_new(node->x_min, mid_y, mid_x, node->y_max);
    node->children[3] = QuadtreeNode_new(mid_x, mid_y, node->x_max, node->y_max);

    for (int i = 0; i < node->line_count; i++) {
      for (int j = 0; j < 4; j++) Quadtree_insert(node->children[j], node->lines[i], depth + 1, dt);
    }
    node->line_count = 0; 
  }
  for (int i = 0; i < 4; i++) Quadtree_insert(node->children[i], l, depth + 1, dt);
}

// O(1) 复杂度的链表合并函数 (用于代替复杂的 Reducer)
void IntersectionEventList_concat(IntersectionEventList* dest, IntersectionEventList* src) {
  if (src->head == NULL) return;
  if (dest->head == NULL) {
    dest->head = src->head;
    dest->tail = src->tail;
  } else {
    dest->tail->next = src->head;
    dest->tail = src->tail;
  }
}

// 【核心并行化优化】使用 cilk_spawn 在四叉树中执行并发深度优先搜索
void Quadtree_getIntersections(QuadtreeNode* node, IntersectionEventList* list, double dt) {
  if (node->is_leaf) {
    for (int i = 0; i < node->line_count; i++) {
      for (int j = i + 1; j < node->line_count; j++) {
        Line *l1 = node->lines[i];
        Line *l2 = node->lines[j];
        if (compareLines(l1, l2) >= 0) {
          Line *temp = l1; l1 = l2; l2 = temp;
        }
        IntersectionType intersectionType = intersect(l1, l2, dt);
        if (intersectionType != NO_INTERSECTION) {
          IntersectionEventList_appendNode(list, l1, l2, intersectionType);
        }
      }
    }
  } else {
    // 为每个并发线程创建独立的局部链表，彻底避免竞态条件 (Race Condition)
    IntersectionEventList local_lists[4];
    for (int i = 0; i < 4; i++) {
      local_lists[i] = IntersectionEventList_make();
    }

    // 衍生子线程并行处理四个象限
    cilk_spawn Quadtree_getIntersections(node->children[0], &local_lists[0], dt);
    cilk_spawn Quadtree_getIntersections(node->children[1], &local_lists[1], dt);
    cilk_spawn Quadtree_getIntersections(node->children[2], &local_lists[2], dt);
    // 最后一个任务直接在当前线程执行，节约调度开销
    Quadtree_getIntersections(node->children[3], &local_lists[3], dt);

    // 等待所有子象限的碰撞检测完成
    cilk_sync;

    // 安全地将局部结果合并回父级链表
    for (int i = 0; i < 4; i++) {
      IntersectionEventList_concat(list, &local_lists[i]);
    }
  }
}

int qsort_compare_events(const void* a, const void* b) {
  IntersectionEventNode* n1 = *(IntersectionEventNode**)a;
  IntersectionEventNode* n2 = *(IntersectionEventNode**)b;
  return IntersectionEventNode_compareData(n1, n2);
}


// =============================================================================
// 原版物理引擎核心逻辑 (含部分并行化优化)
// =============================================================================

CollisionWorld* CollisionWorld_new(const unsigned int capacity) {
  assert(capacity > 0);
  CollisionWorld* collisionWorld = malloc(sizeof(CollisionWorld));
  if (collisionWorld == NULL) return NULL;
  collisionWorld->numLineWallCollisions = 0;
  collisionWorld->numLineLineCollisions = 0;
  collisionWorld->timeStep = 0.5;
  collisionWorld->lines = malloc(capacity * sizeof(Line*));
  collisionWorld->numOfLines = 0;
  return collisionWorld;
}

void CollisionWorld_delete(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) free(collisionWorld->lines[i]);
  free(collisionWorld->lines);
  free(collisionWorld);
}

unsigned int CollisionWorld_getNumOfLines(CollisionWorld* collisionWorld) {
  return collisionWorld->numOfLines;
}

void CollisionWorld_addLine(CollisionWorld* collisionWorld, Line *line) {
  collisionWorld->lines[collisionWorld->numOfLines] = line;
  collisionWorld->numOfLines++;
}

Line* CollisionWorld_getLine(CollisionWorld* collisionWorld, const unsigned int index) {
  if (index >= collisionWorld->numOfLines) return NULL;
  return collisionWorld->lines[index];
}

void CollisionWorld_updateLines(CollisionWorld* collisionWorld) {
  CollisionWorld_detectIntersection(collisionWorld);
  CollisionWorld_updatePosition(collisionWorld);
  CollisionWorld_lineWallCollision(collisionWorld);
}

// 【并行化优化】利用 cilk_for 并发更新每一条线段的位置
void CollisionWorld_updatePosition(CollisionWorld* collisionWorld) {
  double t = collisionWorld->timeStep;
  cilk_for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    line->p1 = Vec_add(line->p1, Vec_multiply(line->velocity, t));
    line->p2 = Vec_add(line->p2, Vec_multiply(line->velocity, t));
  }
}

void CollisionWorld_lineWallCollision(CollisionWorld* collisionWorld) {
  // 墙壁碰撞通常很快，且带有对全局计数器 numLineWallCollisions 的更新
  // 为避免竞态条件，此处保留常规串行 for 循环
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    bool collide = false;

    if ((line->p1.x > BOX_XMAX || line->p2.x > BOX_XMAX) && (line->velocity.x > 0)) {
      line->velocity.x = -line->velocity.x; collide = true;
    }
    if ((line->p1.x < BOX_XMIN || line->p2.x < BOX_XMIN) && (line->velocity.x < 0)) {
      line->velocity.x = -line->velocity.x; collide = true;
    }
    if ((line->p1.y > BOX_YMAX || line->p2.y > BOX_YMAX) && (line->velocity.y > 0)) {
      line->velocity.y = -line->velocity.y; collide = true;
    }
    if ((line->p1.y < BOX_YMIN || line->p2.y < BOX_YMIN) && (line->velocity.y < 0)) {
      line->velocity.y = -line->velocity.y; collide = true;
    }
    if (collide == true) collisionWorld->numLineWallCollisions++;
  }
}

void CollisionWorld_detectIntersection(CollisionWorld* collisionWorld) {
  IntersectionEventList intersectionEventList = IntersectionEventList_make();

  QuadtreeNode* root = QuadtreeNode_new(BOX_XMIN, BOX_YMIN, BOX_XMAX, BOX_YMAX);
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Quadtree_insert(root, collisionWorld->lines[i], 0, collisionWorld->timeStep);
  }

  // 启动并发的深度优先搜索进行碰撞检测
  Quadtree_getIntersections(root, &intersectionEventList, collisionWorld->timeStep);
  QuadtreeNode_delete(root);

  if (intersectionEventList.head != NULL) {
    int count = 0;
    IntersectionEventNode* curr = intersectionEventList.head;
    while (curr != NULL) { count++; curr = curr->next; }

    IntersectionEventNode** array = malloc(count * sizeof(IntersectionEventNode*));
    curr = intersectionEventList.head;
    for (int i = 0; i < count; i++) { array[i] = curr; curr = curr->next; }

    // 使用 qsort 保持处理顺序确定性
    qsort(array, count, sizeof(IntersectionEventNode*), qsort_compare_events);

    intersectionEventList.head = array[0];
    IntersectionEventNode* last = array[0];
    for (int i = 1; i < count; i++) {
      if (IntersectionEventNode_compareData(array[i], last) != 0) {
        last->next = array[i]; last = array[i];
      } else {
        free(array[i]); 
      }
    }
    last->next = NULL;
    free(array);
  }

  // 物理反弹求解必须串行执行，确保物理确定性
  IntersectionEventNode* curNode = intersectionEventList.head;
  while (curNode != NULL) {
    CollisionWorld_collisionSolver(collisionWorld, curNode->l1, curNode->l2, curNode->intersectionType);
    collisionWorld->numLineLineCollisions++;
    curNode = curNode->next;
  }
  IntersectionEventList_deleteNodes(&intersectionEventList);
}

unsigned int CollisionWorld_getNumLineWallCollisions(CollisionWorld* collisionWorld) {
  return collisionWorld->numLineWallCollisions;
}

unsigned int CollisionWorld_getNumLineLineCollisions(CollisionWorld* collisionWorld) {
  return collisionWorld->numLineLineCollisions;
}

void CollisionWorld_collisionSolver(CollisionWorld* collisionWorld, Line *l1, Line *l2, IntersectionType intersectionType) {
  assert(compareLines(l1, l2) < 0);
  assert(intersectionType == L1_WITH_L2 || intersectionType == L2_WITH_L1 || intersectionType == ALREADY_INTERSECTED);

  if (intersectionType == ALREADY_INTERSECTED) {
    Vec p = getIntersectionPoint(l1->p1, l1->p2, l2->p1, l2->p2);
    if (Vec_length(Vec_subtract(l1->p1, p)) < Vec_length(Vec_subtract(l1->p2, p))) {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p2, p)), Vec_length(l1->velocity));
    } else {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p1, p)), Vec_length(l1->velocity));
    }
    if (Vec_length(Vec_subtract(l2->p1, p)) < Vec_length(Vec_subtract(l2->p2, p))) {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p2, p)), Vec_length(l2->velocity));
    } else {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p1, p)), Vec_length(l2->velocity));
    }
    return;
  }

  Vec face, normal;
  if (intersectionType == L1_WITH_L2) {
    Vec v = Vec_makeFromLine(*l2); face = Vec_normalize(v);
  } else {
    Vec v = Vec_makeFromLine(*l1); face = Vec_normalize(v);
  }
  normal = Vec_orthogonal(face);

  double v1Face = Vec_dotProduct(l1->velocity, face), v2Face = Vec_dotProduct(l2->velocity, face);
  double v1Normal = Vec_dotProduct(l1->velocity, normal), v2Normal = Vec_dotProduct(l2->velocity, normal);
  double m1 = Vec_length(Vec_subtract(l1->p1, l1->p2)), m2 = Vec_length(Vec_subtract(l2->p1, l2->p2));
  double newV1Normal = ((m1 - m2) / (m1 + m2)) * v1Normal + (2 * m2 / (m1 + m2)) * v2Normal;
  double newV2Normal = (2 * m1 / (m1 + m2)) * v1Normal + ((m2 - m1) / (m2 + m1)) * v2Normal;

  l1->velocity = Vec_add(Vec_multiply(normal, newV1Normal), Vec_multiply(face, v1Face));
  l2->velocity = Vec_add(Vec_multiply(normal, newV2Normal), Vec_multiply(face, v2Face));
}
