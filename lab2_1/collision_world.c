/** * collision_world.c -- detect and handle line segment intersections
 * Copyright (c) 2012 the Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
 **/

#include "./collision_world.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include "./intersection_detection.h"
#include "./intersection_event_list.h"
#include "./line.h"

// =============================================================================
// 四叉树 (Quadtree) 结构与辅助函数 - 用于优化 $O(N^2)$ 的碰撞检测
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

// 检查线段在时间 dt 内的运动包围盒是否与四叉树节点有重叠
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
  if (!lineOverlapsBox(l, dt, node->x_min, node->y_min, node->x_max, node->y_max)) {
    return;
  }

  if (node->is_leaf) {
    if (node->line_count < node->capacity || depth >= QT_MAX_DEPTH) {
      if (node->line_count == node->capacity) {
        node->capacity *= 2;
        node->lines = realloc(node->lines, sizeof(Line*) * node->capacity);
      }
      node->lines[node->line_count++] = l;
      return;
    }

    // 超过阈值，裂变为 4 个子节点
    node->is_leaf = false;
    double mid_x = (node->x_min + node->x_max) / 2.0;
    double mid_y = (node->y_min + node->y_max) / 2.0;

    node->children[0] = QuadtreeNode_new(node->x_min, node->y_min, mid_x, mid_y);
    node->children[1] = QuadtreeNode_new(mid_x, node->y_min, node->x_max, mid_y);
    node->children[2] = QuadtreeNode_new(node->x_min, mid_y, mid_x, node->y_max);
    node->children[3] = QuadtreeNode_new(mid_x, mid_y, node->x_max, node->y_max);

    for (int i = 0; i < node->line_count; i++) {
      for (int j = 0; j < 4; j++) {
        Quadtree_insert(node->children[j], node->lines[i], depth + 1, dt);
      }
    }
    node->line_count = 0; // 当前节点的线段已全部分配给子节点
  }

  for (int i = 0; i < 4; i++) {
    Quadtree_insert(node->children[i], l, depth + 1, dt);
  }
}

// 递归遍历四叉树，收集所有发生碰撞的对
void Quadtree_getIntersections(QuadtreeNode* node, IntersectionEventList* list, double dt) {
  if (node->is_leaf) {
    for (int i = 0; i < node->line_count; i++) {
      for (int j = i + 1; j < node->line_count; j++) {
        Line *l1 = node->lines[i];
        Line *l2 = node->lines[j];
        if (compareLines(l1, l2) >= 0) {
          Line *temp = l1;
          l1 = l2;
          l2 = temp;
        }
        IntersectionType intersectionType = intersect(l1, l2, dt);
        if (intersectionType != NO_INTERSECTION) {
          IntersectionEventList_appendNode(list, l1, l2, intersectionType);
        }
      }
    }
  } else {
    for (int i = 0; i < 4; i++) {
      Quadtree_getIntersections(node->children[i], list, dt);
    }
  }
}

// 用于快速排序去重的比较函数
int qsort_compare_events(const void* a, const void* b) {
  IntersectionEventNode* n1 = *(IntersectionEventNode**)a;
  IntersectionEventNode* n2 = *(IntersectionEventNode**)b;
  return IntersectionEventNode_compareData(n1, n2);
}


// =============================================================================
// 原版物理引擎核心逻辑 (原封不动保留，确保正确性)
// =============================================================================

CollisionWorld* CollisionWorld_new(const unsigned int capacity) {
  assert(capacity > 0);

  CollisionWorld* collisionWorld = malloc(sizeof(CollisionWorld));
  if (collisionWorld == NULL) {
    return NULL;
  }

  collisionWorld->numLineWallCollisions = 0;
  collisionWorld->numLineLineCollisions = 0;
  collisionWorld->timeStep = 0.5;
  collisionWorld->lines = malloc(capacity * sizeof(Line*));
  collisionWorld->numOfLines = 0;
  return collisionWorld;
}

void CollisionWorld_delete(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    free(collisionWorld->lines[i]);
  }
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

Line* CollisionWorld_getLine(CollisionWorld* collisionWorld,
                             const unsigned int index) {
  if (index >= collisionWorld->numOfLines) {
    return NULL;
  }
  return collisionWorld->lines[index];
}

void CollisionWorld_updateLines(CollisionWorld* collisionWorld) {
  CollisionWorld_detectIntersection(collisionWorld);
  CollisionWorld_updatePosition(collisionWorld);
  CollisionWorld_lineWallCollision(collisionWorld);
}

void CollisionWorld_updatePosition(CollisionWorld* collisionWorld) {
  double t = collisionWorld->timeStep;
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    line->p1 = Vec_add(line->p1, Vec_multiply(line->velocity, t));
    line->p2 = Vec_add(line->p2, Vec_multiply(line->velocity, t));
  }
}

void CollisionWorld_lineWallCollision(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    bool collide = false;

    // Right side
    if ((line->p1.x > BOX_XMAX || line->p2.x > BOX_XMAX)
        && (line->velocity.x > 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Left side
    if ((line->p1.x < BOX_XMIN || line->p2.x < BOX_XMIN)
        && (line->velocity.x < 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Top side
    if ((line->p1.y > BOX_YMAX || line->p2.y > BOX_YMAX)
        && (line->velocity.y > 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Bottom side
    if ((line->p1.y < BOX_YMIN || line->p2.y < BOX_YMIN)
        && (line->velocity.y < 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Update total number of collisions.
    if (collide == true) {
      collisionWorld->numLineWallCollisions++;
    }
  }
}

void CollisionWorld_detectIntersection(CollisionWorld* collisionWorld) {
  IntersectionEventList intersectionEventList = IntersectionEventList_make();

  // 1. 构建四叉树并插入线段
  QuadtreeNode* root = QuadtreeNode_new(BOX_XMIN, BOX_YMIN, BOX_XMAX, BOX_YMAX);
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Quadtree_insert(root, collisionWorld->lines[i], 0, collisionWorld->timeStep);
  }

  // 2. 使用四叉树收集碰撞对
  Quadtree_getIntersections(root, &intersectionEventList, collisionWorld->timeStep);
  QuadtreeNode_delete(root);

  // 3. 将链表转为数组，使用 qsort 极速排序并去重
  if (intersectionEventList.head != NULL) {
    int count = 0;
    IntersectionEventNode* curr = intersectionEventList.head;
    while (curr != NULL) {
      count++;
      curr = curr->next;
    }

    IntersectionEventNode** array = malloc(count * sizeof(IntersectionEventNode*));
    curr = intersectionEventList.head;
    for (int i = 0; i < count; i++) {
      array[i] = curr;
      curr = curr->next;
    }

    // 采用与原版选择排序完全相同的比较规则进行快排
    qsort(array, count, sizeof(IntersectionEventNode*), qsort_compare_events);

    // 去重逻辑：跨界线段可能会被插入多个象限，产生相同的重复碰撞事件
    intersectionEventList.head = array[0];
    IntersectionEventNode* last = array[0];
    for (int i = 1; i < count; i++) {
      if (IntersectionEventNode_compareData(array[i], last) != 0) {
        last->next = array[i];
        last = array[i];
      } else {
        free(array[i]); // 释放重复项，防止计数错误和二次反弹
      }
    }
    last->next = NULL;
    free(array);
  }

  // 4. 调用原版碰撞处理器并计数
  IntersectionEventNode* curNode = intersectionEventList.head;
  while (curNode != NULL) {
    CollisionWorld_collisionSolver(collisionWorld, curNode->l1, curNode->l2,
                                   curNode->intersectionType);
    collisionWorld->numLineLineCollisions++;
    curNode = curNode->next;
  }

  IntersectionEventList_deleteNodes(&intersectionEventList);
}

unsigned int CollisionWorld_getNumLineWallCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineWallCollisions;
}

unsigned int CollisionWorld_getNumLineLineCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineLineCollisions;
}

void CollisionWorld_collisionSolver(CollisionWorld* collisionWorld,
                                    Line *l1, Line *l2,
                                    IntersectionType intersectionType) {
  assert(compareLines(l1, l2) < 0);
  assert(intersectionType == L1_WITH_L2
         || intersectionType == L2_WITH_L1
         || intersectionType == ALREADY_INTERSECTED);

  // Despite our efforts to determine whether lines will intersect ahead
  // of time (and to modify their velocities appropriately), our
  // simplified model can sometimes cause lines to intersect.  In such a
  // case, we compute velocities so that the two lines can get unstuck in
  // the fastest possible way, while still conserving momentum and kinetic
  // energy.
  if (intersectionType == ALREADY_INTERSECTED) {
    Vec p = getIntersectionPoint(l1->p1, l1->p2, l2->p1, l2->p2);

    if (Vec_length(Vec_subtract(l1->p1, p))
        < Vec_length(Vec_subtract(l1->p2, p))) {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p2, p)),
                                  Vec_length(l1->velocity));
    } else {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p1, p)),
                                  Vec_length(l1->velocity));
    }
    if (Vec_length(Vec_subtract(l2->p1, p))
        < Vec_length(Vec_subtract(l2->p2, p))) {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p2, p)),
                                  Vec_length(l2->velocity));
    } else {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p1, p)),
                                  Vec_length(l2->velocity));
    }
    return;
  }

  // Compute the collision face/normal vectors.
  Vec face;
  Vec normal;
  if (intersectionType == L1_WITH_L2) {
    Vec v = Vec_makeFromLine(*l2);
    face = Vec_normalize(v);
  } else {
    Vec v = Vec_makeFromLine(*l1);
    face = Vec_normalize(v);
  }
  normal = Vec_orthogonal(face);

  // Obtain each line's velocity components with respect to the collision
  // face/normal vectors.
  double v1Face = Vec_dotProduct(l1->velocity, face);
  double v2Face = Vec_dotProduct(l2->velocity, face);
  double v1Normal = Vec_dotProduct(l1->velocity, normal);
  double v2Normal = Vec_dotProduct(l2->velocity, normal);

  // Compute the mass of each line (we simply use its length).
  double m1 = Vec_length(Vec_subtract(l1->p1, l1->p2));
  double m2 = Vec_length(Vec_subtract(l2->p1, l2->p2));

  // Perform the collision calculation (computes the new velocities along
  // the direction normal to the collision face such that momentum and
  // kinetic energy are conserved).
  double newV1Normal = ((m1 - m2) / (m1 + m2)) * v1Normal
      + (2 * m2 / (m1 + m2)) * v2Normal;
  double newV2Normal = (2 * m1 / (m1 + m2)) * v1Normal
      + ((m2 - m1) / (m2 + m1)) * v2Normal;

  // Combine the resulting velocities.
  l1->velocity = Vec_add(Vec_multiply(normal, newV1Normal),
                         Vec_multiply(face, v1Face));
  l2->velocity = Vec_add(Vec_multiply(normal, newV2Normal),
                         Vec_multiply(face, v2Face));

  return;
}
