#pragma once


//#pragma warning (disable:4311)

#include <vector>
#include <list>
#include "point2d.h"
#include "point3d.h"
#include "point4d.h"

typedef unsigned int UInt;

#ifndef NOMINMAX
#define NOMINMAX
#endif
#undef min
//#define STORE_CELLS_IN_TREE

namespace JBSlib {

// ****************************************************************************************

template <class T> class DFOcTree;

// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

//! Cell of a n-D Octree.
template <class T>
class DFOcNode {

public:

// ****************************************************************************************

	typedef T Point;
	typedef DFOcNode<Point> OcNodeType;
	typedef typename T::ValueType ValueType;
	static const UInt dim = Point::dim;
	//! Number of children = 2_tothepowerof_dim (8 for 3D, 16 for 4D...).
	static const UInt numChildren = 1 << dim;

// ****************************************************************************************

	DFOcNode(DFOcTree< OcNodeType >* tree_, Point min_, Point max_, OcNodeType* father_ = NULL): min_ext(min_), max_ext(max_), tree(tree_), m_isLeafNode(true), father(father_) {
		if (father == NULL) nodelevel = 0;
		else nodelevel = father->getNodeLevel()+1;

		center = (min_ext+max_ext)/2;
	}

// ****************************************************************************************

	~DFOcNode() {
		if (!m_isLeafNode) { // we have kids, kill them!
			for (UInt i1 = 0; i1 < numChildren; i1++) {
				delete children[i1];
			}
			delete[] children;
		}
	}

// ****************************************************************************************

	//! Prints the content of the current node.
	inline void printContent() const {
		std::cerr << "Node: " << std::endl;
		for (UInt i1 = 0; i1 < content.size(); i1++) {
			//content[i1]->print();
			std::cerr << content[i1] << " ";
		}
		std::cerr << std::endl;
		std::cerr << "------------------" << std::endl;
	}

// ****************************************************************************************

	//! Returns true if the node does not contain any further kids.
	inline bool isLeafNode() const {
		return m_isLeafNode;
	}
	
// ****************************************************************************************

	//! Returns the number of triangles in this node.
	inline UInt getNumObjects() const {
		return (UInt)content.size();
	}


// ****************************************************************************************

	//! Returns the min corner of the cell.
	Point getMin() const {
		return min_ext;
	}

// ****************************************************************************************

	//! Returns the min corner of the cell.
	Point getMax() const {
		return max_ext;
	}

// ****************************************************************************************

	UInt getNodeLevel() const {
		return nodelevel;
	}

// ****************************************************************************************

	//! Inserts an object into this node.
	inline int insertObject(UInt obj_pos) {

		int pos = 0;

		if (m_isLeafNode) {

			// Check if the object is already contained in the current node
			for (UInt i1 = 0; i1 < content.size(); i1++) {
				//if (tree->objects[obj_pos] == tree->objects[content[i1]]) {
				//	return content[i1];
				//}
				if ((tree->objects[obj_pos].dist(tree->objects[content[i1]]) < std::numeric_limits<float>::min())) {
					return content[i1];
				}
			}

			if (getNumObjects() >= tree->getMaxNumObjectsPerNode()) {

				// split node;
				splitNode();

				// distribute old objects to kids:
				for (UInt i1 = 0; i1 < content.size(); i1++) {
					propagateToChildrenAvoidMultiple(content[i1]);
				}
				// Add newly added point to children.
				pos = propagateToChildrenAvoidMultiple(obj_pos);

				// remove own storage:
				content.clear();

			} else { // We have still some slots for points left.
				content.push_back(obj_pos);
				return -1;
			}
		} else { // We're not a leaf-node, let out kids handle the point
			// propagate object to children.
			pos = propagateToChildrenAvoidMultiple(obj_pos);
		}
		return pos;
	}

// ****************************************************************************************

protected:

// ****************************************************************************************

	inline Point getInterpolatedPos(UInt x, UInt y) const {
		Point p;
		p[0] = (min_ext.x*(2-x) + max_ext.x*(x)) / ValueType(2.0);
		p[1] = (min_ext.y*(2-y) + max_ext.y*(y)) / ValueType(2.0);
		return p;
	}

// ****************************************************************************************

	inline Point getInterpolatedPos(UInt x, UInt y, UInt z) const {
		Point p;
		p[0] = (min_ext.x*(2-x) + max_ext.x*(x)) / ValueType(2.0);
		p[1] = (min_ext.y*(2-y) + max_ext.y*(y)) / ValueType(2.0);
		p[2] = (min_ext.z*(2-z) + max_ext.z*(z)) / ValueType(2.0);
		return p;
	}

// ****************************************************************************************

	inline Point getInterpolatedPos(UInt x, UInt y, UInt z, UInt w) const {
		Point p;
		p[0] = (min_ext.x*(2-x) + max_ext.x*(x)) / ValueType(2.0);
		p[1] = (min_ext.y*(2-y) + max_ext.y*(y)) / ValueType(2.0);
		p[2] = (min_ext.z*(2-z) + max_ext.z*(z)) / ValueType(2.0);
		p[3] = (min_ext.z*(2-w) + max_ext.z*(w)) / ValueType(2.0);
		return p;
	}

// ****************************************************************************************

	int propagateToChildrenAvoidMultiple(UInt obj_pos);
	
// ****************************************************************************************

	int propagateToChildrenAvoidMultiple2D(UInt obj_pos) {
		const Point& obj = tree->objects[obj_pos];
		int xx = (obj[0] > center[0]); int yy = (obj[1] > center[1]);
		int pos = children[getPosInChildrenArray(xx,yy)]->insertObject(obj_pos);
		return pos;
	}

// ****************************************************************************************

	int propagateToChildrenAvoidMultiple3D(UInt obj_pos) {
		const Point& obj = tree->objects[obj_pos];
		int xx = (obj[0] > center[0]); int yy = (obj[1] > center[1]); int zz = (obj[2] > center[2]);
		int pos = children[getPosInChildrenArray(xx,yy,zz)]->insertObject(obj_pos);
		return pos;
	}

// ****************************************************************************************

	int propagateToChildrenAvoidMultiple4D(UInt obj_pos) {
		const Point& obj = tree->objects[obj_pos];
		int xx = (obj[0] > center[0]); int yy = (obj[1] > center[1]); int zz = (obj[2] > center[2]); int ww = (obj[3] > center[3]);
		int pos = children[getPosInChildrenArray(xx,yy,zz,ww)]->insertObject(obj_pos);
		return pos;
	}

// ****************************************************************************************

	void splitNode();

// ****************************************************************************************

	void splitNode2D() {
		children = new DFOcNode<Point>*[numChildren];
		children[getPosInChildrenArray(0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0), getInterpolatedPos(1,1), this);
		children[getPosInChildrenArray(1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0), getInterpolatedPos(2,1), this);
		children[getPosInChildrenArray(0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1), getInterpolatedPos(1,2), this);
		children[getPosInChildrenArray(1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1), getInterpolatedPos(2,2), this);
		// We're not longer leaf-node:
		m_isLeafNode = false;
	}

// ****************************************************************************************

	void splitNode3D() {
		children = new DFOcNode<Point>*[numChildren];
		children[getPosInChildrenArray(0,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,0), getInterpolatedPos(1,1,1), this);
		children[getPosInChildrenArray(1,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,0), getInterpolatedPos(2,1,1), this);
		children[getPosInChildrenArray(0,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,0), getInterpolatedPos(1,2,1), this);
		children[getPosInChildrenArray(1,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,0), getInterpolatedPos(2,2,1), this);
		children[getPosInChildrenArray(0,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,1), getInterpolatedPos(1,1,2), this);
		children[getPosInChildrenArray(1,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,1), getInterpolatedPos(2,1,2), this);
		children[getPosInChildrenArray(0,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,1), getInterpolatedPos(1,2,2), this);
		children[getPosInChildrenArray(1,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,1), getInterpolatedPos(2,2,2), this);
		// We're not longer leaf-node:
		m_isLeafNode = false;
	}

// ****************************************************************************************

	void splitNode4D() {
		children = new DFOcNode<Point>*[numChildren];
		children[getPosInChildrenArray(0,0,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,0,0), getInterpolatedPos(1,1,1,1), this);
		children[getPosInChildrenArray(1,0,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,0,0), getInterpolatedPos(2,1,1,1), this);
		children[getPosInChildrenArray(0,1,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,0,0), getInterpolatedPos(1,2,1,1), this);
		children[getPosInChildrenArray(1,1,0,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,0,0), getInterpolatedPos(2,2,1,1), this);
		children[getPosInChildrenArray(0,0,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,1,0), getInterpolatedPos(1,1,2,1), this);
		children[getPosInChildrenArray(1,0,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,1,0), getInterpolatedPos(2,1,2,1), this);
		children[getPosInChildrenArray(0,1,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,1,0), getInterpolatedPos(1,2,2,1), this);
		children[getPosInChildrenArray(1,1,1,0)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,1,0), getInterpolatedPos(2,2,2,1), this);
		children[getPosInChildrenArray(0,0,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,0,1), getInterpolatedPos(1,1,1,2), this);
		children[getPosInChildrenArray(1,0,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,0,1), getInterpolatedPos(2,1,1,2), this);
		children[getPosInChildrenArray(0,1,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,0,1), getInterpolatedPos(1,2,1,2), this);
		children[getPosInChildrenArray(1,1,0,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,0,1), getInterpolatedPos(2,2,1,2), this);
		children[getPosInChildrenArray(0,0,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,0,1,1), getInterpolatedPos(1,1,2,2), this);
		children[getPosInChildrenArray(1,0,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,0,1,1), getInterpolatedPos(2,1,2,2), this);
		children[getPosInChildrenArray(0,1,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(0,1,1,1), getInterpolatedPos(1,2,2,2), this);
		children[getPosInChildrenArray(1,1,1,1)] = new DFOcNode<T>(tree, getInterpolatedPos(1,1,1,1), getInterpolatedPos(2,2,2,2), this);
		// We're not longer leaf-node:
		m_isLeafNode = false;
	}

// ****************************************************************************************

	inline UInt getPosInChildrenArray(UInt x, UInt y) const					{ return x*2+y; }
	inline UInt getPosInChildrenArray(UInt x, UInt y, UInt z) const			{ return x*4+y*2+z; }
	inline UInt getPosInChildrenArray(UInt x, UInt y, UInt z, UInt w) const	{ return x*8+y*4+z*2+w; }

// ****************************************************************************************

	//! Depth of this node
	UInt nodelevel;

	//! Min. corner of cell.
	Point min_ext;
	//! Max. corner of cell.
	Point max_ext;
	//! Center of the cell = (max_ext+min_ext)/2;
	Point center;

	//! 'Father'-tree.
	DFOcTree< DFOcNode<Point> >* tree;

	//! Flag for leaf-/inner-nodes.
	bool m_isLeafNode;

	//! Children nodes.
	OcNodeType** children;
	//! Father-node pointer.
	OcNodeType* father;

	//! Indices (in parent-tree array) of the cells content.
	std::vector<UInt> content;

};

// ****************************************************************************************
// ************************** TEMPLATE SPECIFIC FUNCTIONS *********************************
// ****************************************************************************************

template <> void DFOcNode<vec2d>::splitNode() { splitNode2D(); }
template <> void DFOcNode<vec2f>::splitNode() { splitNode2D(); }
template <> void DFOcNode<vec3d>::splitNode() { splitNode3D(); }
template <> void DFOcNode<vec3f>::splitNode() { splitNode3D(); }
template <> void DFOcNode<vec4d>::splitNode() { splitNode4D(); }
template <> void DFOcNode<vec4f>::splitNode() { splitNode4D(); }
template <> int DFOcNode<vec2d>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple2D(obj_pos); }
template <> int DFOcNode<vec2f>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple2D(obj_pos); }
template <> int DFOcNode<vec3d>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple3D(obj_pos); }
template <> int DFOcNode<vec3f>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple3D(obj_pos); }
template <> int DFOcNode<vec4d>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple4D(obj_pos); }
template <> int DFOcNode<vec4f>::propagateToChildrenAvoidMultiple(UInt obj_pos) { return propagateToChildrenAvoidMultiple4D(obj_pos); }



typedef DFOcNode<vec2f> DFOcNode2f;
typedef DFOcNode<vec2d> DFOcNode2d;
typedef DFOcNode<vec3f> DFOcNode3f;
typedef DFOcNode<vec3d> DFOcNode3d;
typedef DFOcNode<vec4f> DFOcNode4f;
typedef DFOcNode<vec4d> DFOcNode4d;



// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************







//! An nD-Octree for detecting multiple points.
/*!
	use DFOcNode<point2d>, DFOcNode<point3d> or DFOcNode<point4d> as template parameter.
*/
template <class T>
class DFOcTree {

public:

// ****************************************************************************************

	typedef T OcNodeType;
	typedef typename T::Point Point;
	typedef typename Point::ValueType ValueType;
	static const UInt dim = Point::dim;

	//! Allow DFOcNode to access our private members.
	friend class DFOcNode<Point>;

// ****************************************************************************************

	//! Allocate a new octree ranging from 'min_' to 'max_'.
	DFOcTree(UInt maxNumObjectsPerNode_ = 10, Point min_ = Point((ValueType)0), Point max_ = Point((ValueType)1))
	: min_ext(min_), max_ext(max_), maxNumObjectsPerNode(maxNumObjectsPerNode_)
	{
		//objects_vector.reserve(8000000);
		head = new OcNodeType(this, min_ext, max_ext);
		count = 0;
		pointsDuplicate = 0;
		pointsAdded = 0;

		//std::cerr << "Created octree for: (" << min_ext << ") to (" << max_ext << ")" << std::endl;

	}

// ****************************************************************************************

	//! Destructor.
	~DFOcTree() {
		delete head;
	}

// ****************************************************************************************

	//! Inserts an object into the octree. returns false if the point was already contained in the octree.
	inline bool insertObject(UInt& diff, const Point& obj_) {

// Define 'PERFORM_NAN_TESTS' in your main file before including this header if you wish to check for 'nan's!
#ifdef PERFORM_NAN_TESTS
		for (UInt i = 0; i < Point::dim; ++i) {
			if (obj_[i] != obj_[i]) {
				std::cerr << "Not a number!" << std::endl;
				diff = 0;
				return false;
			}
		}
#endif
		
		pointsAdded++;

		objects.push_back(obj_);

		int pos = head->insertObject(count);

		if (pos < 0) { // Point is unique.
			diff = count;
			count++;
			return true;
		}

		//! Remove duplicate again.
		objects.pop_back();

		pointsDuplicate++;

		diff = pos;

		return false;
	}

// ****************************************************************************************

	//! Returns how many objects may be maximally stored in each leaf-node.
	inline UInt getMaxNumObjectsPerNode() const {
		return maxNumObjectsPerNode;
	}

// ****************************************************************************************

	//! Returns the 'min' corner of the octree.
	inline Point getMin() const {
		return min_ext;
	}

// ****************************************************************************************

	//! Returns the 'max' corner of the octree.
	inline Point getMax() const {
		return max_ext;
	}

// ****************************************************************************************

	//! Number of (unique) points contained in tree.
	UInt count;

	//! Number of points that have not been added because they were duplicates.
	UInt pointsDuplicate;

	//! Total number of points that have been added (= 'count' + 'pointsDuplicate').
	UInt pointsAdded;

	//! Returns an 'Point*' to the points in the octree
	Point* getObjects() {
		return &objects[0];
	}

protected:

// ****************************************************************************************

	//! Storage for points in tree.
	std::vector<Point> objects;

	//! Min. corner of head-cell.
	Point min_ext;
	//! Max. corner of head-cell.
	Point max_ext;

	//! Max. number of objects in each cell before cell is splitted.
	UInt maxNumObjectsPerNode;

	//! First cell.
	OcNodeType* head;
};

typedef DFOcTree< DFOcNode2f > DFOctree2f;
typedef DFOcTree< DFOcNode2d > DFOctree2d;
typedef DFOcTree< DFOcNode3f > DFOctree3f;
typedef DFOcTree< DFOcNode3d > DFOctree3d;
typedef DFOcTree< DFOcNode4f > DFOctree4f;
typedef DFOcTree< DFOcNode4d > DFOctree4d;

}