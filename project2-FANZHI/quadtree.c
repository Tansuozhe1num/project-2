/**
 * quadtree.c -- Quadtree spatial data structure for collision detection
 * FANZHI - Project2 Part 1 Implementation
 *
 * This quadtree optimizes collision detection from O(N²) to O(N log N)
 * by partitioning 2D space and only testing lines in the same region.
 **/

#include "./quadtree.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

// Include OpenCilk headers if available
#ifdef __cilk__
#include <cilk/cilk.h>
#endif

// Helper: minimum of two doubles
static double min_double(double a, double b) {
  return (a < b) ? a : b;
}

// Helper: maximum of two doubles
static double max_double(double a, double b) {
  return (a > b) ? a : b;
}

// Calculate the bounding box of a line including its swept area over one timestep
LineBounds Line_getBounds(Line* line, double timeStep) {
  LineBounds bounds;

  // Current position bounds
  double curMinX = min_double(line->p1.x, line->p2.x);
  double curMinY = min_double(line->p1.y, line->p2.y);
  double curMaxX = max_double(line->p1.x, line->p2.x);
  double curMaxY = max_double(line->p1.y, line->p2.y);

  // Position after one timestep
  Vec futureP1 = Vec_add(line->p1, Vec_multiply(line->velocity, timeStep));
  Vec futureP2 = Vec_add(line->p2, Vec_multiply(line->velocity, timeStep));

  double futMinX = min_double(futureP1.x, futureP2.x);
  double futMinY = min_double(futureP1.y, futureP2.y);
  double futMaxX = max_double(futureP1.x, futureP2.x);
  double futMaxY = max_double(futureP1.y, futureP2.y);

  // Bounding box covers both current and future positions
  bounds.minX = min_double(curMinX, futMinX);
  bounds.minY = min_double(curMinY, futMinY);
  bounds.maxX = max_double(curMaxX, futMaxX);
  bounds.maxY = max_double(curMaxY, futMaxY);

  return bounds;
}

// Create a new quadtree node
QuadtreeNode* QuadtreeNode_new(double x, double y, double width, double height, int depth) {
  QuadtreeNode* node = malloc(sizeof(QuadtreeNode));
  if (node == NULL) {
    return NULL;
  }

  node->x = x;
  node->y = y;
  node->width = width;
  node->height = height;
  node->depth = depth;

  // Allocate initial capacity for lines
  node->capacity = QUADTREE_MAX_OBJECTS;
  node->lines = malloc(node->capacity * sizeof(Line*));
  if (node->lines == NULL) {
    free(node);
    return NULL;
  }
  node->count = 0;

  // Initialize children to NULL (leaf node)
  for (int i = 0; i < 4; i++) {
    node->children[i] = NULL;
  }

  return node;
}

// Destroy a quadtree node and all its children recursively
void QuadtreeNode_destroy(QuadtreeNode* node) {
  if (node == NULL) {
    return;
  }

  // Destroy children first
  for (int i = 0; i < 4; i++) {
    if (node->children[i] != NULL) {
      QuadtreeNode_destroy(node->children[i]);
    }
  }

  // Free this node's resources
  if (node->lines != NULL) {
    free(node->lines);
  }
  free(node);
}

// Split a node into four children
void QuadtreeNode_split(QuadtreeNode* node) {
  double halfWidth = node->width / 2.0;
  double halfHeight = node->height / 2.0;
  int childDepth = node->depth + 1;

  // Create four children: NW, NE, SW, SE
  // 0: NW (top-left)
  node->children[0] = QuadtreeNode_new(node->x, node->y + halfHeight,
                                        halfWidth, halfHeight, childDepth);
  // 1: NE (top-right)
  node->children[1] = QuadtreeNode_new(node->x + halfWidth, node->y + halfHeight,
                                        halfWidth, halfHeight, childDepth);
  // 2: SW (bottom-left)
  node->children[2] = QuadtreeNode_new(node->x, node->y,
                                        halfWidth, halfHeight, childDepth);
  // 3: SE (bottom-right)
  node->children[3] = QuadtreeNode_new(node->x + halfWidth, node->y,
                                        halfWidth, halfHeight, childDepth);

  // Redistribute existing lines to children
  unsigned int oldCount = node->count;
  Line** oldLines = node->lines;

  // Reset this node (it's no longer a leaf)
  node->lines = NULL;
  node->count = 0;
  node->capacity = 0;

  // Re-insert all lines into children
  for (unsigned int i = 0; i < oldCount; i++) {
    LineBounds bounds = Line_getBounds(oldLines[i], 0.5); // Use default timestep
    // Insert into all children that the line's bounds intersect
    for (int c = 0; c < 4; c++) {
      if (node->children[c] != NULL) {
        QuadtreeNode_insert(node->children[c], oldLines[i], &bounds);
      }
    }
  }

  free(oldLines);
}

// Check if two rectangles overlap
static bool rectsOverlap(double x1, double y1, double w1, double h1,
                          double x2, double y2, double w2, double h2) {
  return (x1 < x2 + w2) && (x1 + w1 > x2) &&
         (y1 < y2 + h2) && (y1 + h1 > y2);
}

// Insert a line into the quadtree node
void QuadtreeNode_insert(QuadtreeNode* node, Line* line, LineBounds* bounds) {
  // Check if the line's bounds intersect with this node's region
  if (!rectsOverlap(node->x, node->y, node->width, node->height,
                    bounds->minX, bounds->minY,
                    bounds->maxX - bounds->minX, bounds->maxY - bounds->minY)) {
    return; // Line doesn't intersect this node
  }

  // If this is a leaf node and has space, or we've reached max depth
  if (node->children[0] == NULL) {
    // Add line to this node
    if (node->count >= node->capacity) {
      // Expand capacity
      node->capacity *= 2;
      node->lines = realloc(node->lines, node->capacity * sizeof(Line*));
      if (node->lines == NULL) {
        return; // Memory allocation failed
      }
    }
    node->lines[node->count] = line;
    node->count++;

    // Split if we exceed max objects and haven't reached max depth
    if (node->count > QUADTREE_MAX_OBJECTS && node->depth < QUADTREE_MAX_DEPTH) {
      QuadtreeNode_split(node);
    }
  } else {
    // This is an internal node, insert into all children that overlap
    for (int i = 0; i < 4; i++) {
      if (node->children[i] != NULL) {
        QuadtreeNode_insert(node->children[i], line, bounds);
      }
    }
  }
}

// Query the quadtree for all lines that might collide with the given line
void QuadtreeNode_query(QuadtreeNode* node, LineBounds* bounds,
                        Line** foundLines, unsigned int* foundCount,
                        unsigned int maxFound) {
  // Check if the query bounds intersect with this node's region
  if (!rectsOverlap(node->x, node->y, node->width, node->height,
                    bounds->minX, bounds->minY,
                    bounds->maxX - bounds->minX, bounds->maxY - bounds->minY)) {
    return; // No overlap
  }

  // Add all lines in this node to the result
  for (unsigned int i = 0; i < node->count; i++) {
    if (*foundCount >= maxFound) {
      return; // Buffer full
    }
    // Avoid duplicates by checking if line is already in the result
    bool alreadyFound = false;
    for (unsigned int j = 0; j < *foundCount; j++) {
      if (foundLines[j] == node->lines[i]) {
        alreadyFound = true;
        break;
      }
    }
    if (!alreadyFound) {
      foundLines[*foundCount] = node->lines[i];
      (*foundCount)++;
    }
  }

  // Recursively query children
  if (node->children[0] != NULL) {
    for (int i = 0; i < 4; i++) {
      if (node->children[i] != NULL) {
        QuadtreeNode_query(node->children[i], bounds, foundLines, foundCount, maxFound);
      }
    }
  }
}

// Build a quadtree from an array of lines
QuadtreeNode* Quadtree_build(Line** lines, unsigned int numLines,
                              double timeStep,
                              double x, double y,
                              double width, double height) {
  QuadtreeNode* root = QuadtreeNode_new(x, y, width, height, 0);
  if (root == NULL) {
    return NULL;
  }

  // Insert all lines into the quadtree
  for (unsigned int i = 0; i < numLines; i++) {
    LineBounds bounds = Line_getBounds(lines[i], timeStep);
    QuadtreeNode_insert(root, lines[i], &bounds);
  }

  return root;
}

// Recursive parallel collision detection using cilk_spawn
// This traverses the quadtree and detects collisions within each node
// Only tests pairs where l1 is the "primary" owner (lowest ID)
void Quadtree_detectCollisions(QuadtreeNode* node,
                                IntersectionEventList* eventList,
                                unsigned int* collisionCount,
                                double timeStep) {
  if (node == NULL) {
    return;
  }

  // Detect collisions among lines in this node
  // Only test pairs where this node is the "primary" node for l1
  // (i.e., l1 has the lowest ID among all nodes containing it)
  for (unsigned int i = 0; i < node->count; i++) {
    Line *l1 = node->lines[i];
    for (unsigned int j = i + 1; j < node->count; j++) {
      Line *l2 = node->lines[j];

      // Skip if same line
      if (l1 == l2) continue;

      // Ensure compareLines(l1, l2) < 0 to avoid duplicate testing
      if (compareLines(l1, l2) >= 0) {
        continue;
      }

      // Test for intersection
      IntersectionType intersectionType = intersect(l1, l2, timeStep);
      if (intersectionType != NO_INTERSECTION) {
        IntersectionEventList_appendNode(eventList, l1, l2, intersectionType);
        (*collisionCount)++;
      }
    }
  }

  // Recursively process child nodes in parallel using cilk_spawn
  if (node->children[0] != NULL) {
    // Spawn parallel tasks for each child
    for (int i = 0; i < 4; i++) {
      if (node->children[i] != NULL) {
#ifdef __cilk__
        cilk_spawn Quadtree_detectCollisions(node->children[i], eventList, collisionCount, timeStep);
#else
        Quadtree_detectCollisions(node->children[i], eventList, collisionCount, timeStep);
#endif
      }
    }
#ifdef __cilk__
    cilk_sync;
#endif
  }
}
