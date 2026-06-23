/**
 * quadtree.h -- Quadtree spatial data structure for collision detection
 * FANZHI - Project2 Part 1 Implementation
 *
 * This quadtree optimizes collision detection from O(N²) to O(N log N)
 * by partitioning 2D space and only testing lines in the same region.
 **/

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include "./line.h"
#include "./intersection_event_list.h"

// Quadtree parameters
#define QUADTREE_MAX_OBJECTS 10
#define QUADTREE_MAX_DEPTH 8

// Bounding box for a line (includes swept area over one timestep)
struct LineBounds {
  double minX;
  double minY;
  double maxX;
  double maxY;
};
typedef struct LineBounds LineBounds;

// Quadtree node structure
struct QuadtreeNode {
  double x;           // Left edge of this node's region
  double y;           // Bottom edge of this node's region
  double width;       // Width of this node's region
  double height;      // Height of this node's region

  Line** lines;       // Array of line pointers stored in this node
  unsigned int count;         // Number of lines currently stored
  unsigned int capacity;      // Capacity of the lines array

  int depth;          // Current depth in the tree (root = 0)

  struct QuadtreeNode* children[4]; // Four child nodes (NULL if leaf)
  // Children indices: 0=NW, 1=NE, 2=SW, 3=SE
};
typedef struct QuadtreeNode QuadtreeNode;

// Calculate the bounding box of a line including its swept area over one timestep
LineBounds Line_getBounds(Line* line, double timeStep);

// Create a new quadtree node covering the given region
QuadtreeNode* QuadtreeNode_new(double x, double y, double width, double height, int depth);

// Destroy a quadtree node and all its children recursively
void QuadtreeNode_destroy(QuadtreeNode* node);

// Split a node into four children
void QuadtreeNode_split(QuadtreeNode* node);

// Determine which child quadrant(s) a line's bounds intersect with
// Returns a bitmask: bit 0=NW, bit 1=NE, bit 2=SW, bit 3=SE
int QuadtreeNode_getChildIndex(QuadtreeNode* node, LineBounds* bounds);

// Insert a line into the quadtree node
// The line may be inserted into multiple children if it spans quadrants
void QuadtreeNode_insert(QuadtreeNode* node, Line* line, LineBounds* bounds);

// Query the quadtree for all lines that might collide with the given line
// Returns lines from the same node and ancestor nodes
// foundLines: output array to store found lines
// foundCount: output count of found lines
// maxFound: maximum number of lines to find
void QuadtreeNode_query(QuadtreeNode* node, LineBounds* bounds,
                        Line** foundLines, unsigned int* foundCount,
                        unsigned int maxFound);

// Build a quadtree from an array of lines
QuadtreeNode* Quadtree_build(Line** lines, unsigned int numLines,
                              double timeStep,
                              double x, double y,
                              double width, double height);

// Recursive parallel collision detection using cilk_spawn
// This traverses the quadtree and detects collisions within each node
void Quadtree_detectCollisions(QuadtreeNode* node,
                                IntersectionEventList* eventList,
                                unsigned int* collisionCount,
                                double timeStep);

#endif // QUADTREE_H_
