#ifndef MESH_H_
#define MESH_H_

#include "fdaPDE.h"
#include "mesh_objects.h"



using std::vector;

template <UInt ORDER,UInt mydim, UInt ndim>
class MeshHandler{
};

//!  2D MESH:
//!  This class gives an object-oriented reading interface to the output of the library Triangle (Jonathan Richard Shewchuk).
/*!
 * The template parameters specify the order of its elements.
 * The aim of this class is to do not introduce any initialization overhead,
 * beacuse it will be called many time during the execution of a R script
*/
template <UInt ORDER>
class MeshHandler<ORDER,2,2> {
public:
	typedef int UInt;
	//! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
      * constructed with the TriLibrary (our R wrapper for the Triangle library)
    */

    MeshHandler(Real* points, UInt* edges, UInt* triangles, UInt* neighbors, UInt num_nodes, UInt num_edges, UInt num_triangles):
			points_(points), edges_(edges), triangles_(triangles), neighbors_(neighbors), num_nodes_(num_nodes), num_edges_(num_edges), num_triangles_(num_triangles) {};

    #ifdef R_VERSION_
	MeshHandler(SEXP Rmesh);
	#endif

	~MeshHandler(){};

	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_nodes() const {return num_nodes_;}

	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_triangles() const {return num_triangles_;}

    //! A normal member returning an unsigned integer value.
    /*!
      \return The number of edges in the mesh
    */
    UInt num_edges() const {return num_edges_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Edge
    /*!
     * \param id an Id argument
      \return The edge with the specified id
    */
    Edge getEdge(Id id);

   //! A normal member setting a Triangle
        /*!
         * \param id an Id argument
          \return The triangle with order coerent to that of the mesh with the specified id
        */
    //void setTriangle(Triangle<ORDER * 3,2,2>& tri, Id id) const;

    //! A normal member returning a Triangle
    /*!
     * \param id an Id argument 
      \return The triangle with order coerent to that of the mesh with the specified id
    */ 
    Triangle<ORDER * 3,2,2>  getTriangle(Id id) const;

    //The "number" neighbor of triangle i is opposite the "number" corner of triangle i
    //! A normal member returning the Neighbors of a triangle
    /*!
     * \param id the id of the triangle
     * \param number the number of the vertex
      \return The triangle that has as an edge the one opposite to the specified
      vertex
    */
    Triangle<ORDER * 3,2,2> getNeighbors(Id id_triangle, UInt number) const;

    void printPoints(std::ostream & out);
    void printEdges(std::ostream & out);
    void printTriangles(std::ostream & out);
    void printNeighbors(std::ostream & out);

     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a simply research between all triangle of the mesh
     * \param point the point we want to locate
      \return The triangle that contains the point
    */
    Triangle<ORDER * 3,2,2> findLocationNaive(Point point) const;

     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
     * \param point the point we want to locate
     * \param starting_triangles a vector of points that specifies the poposed starting
     * points for the walking algorithm
      \return The triangle that contains the point
    */
    Triangle<ORDER * 3,2,2> findLocationWalking(const Point& point, const Triangle<ORDER * 3,2,2>& starting_triangle) const;

    //int readMesh(std::string const & file);
	//double measure()const;
	//bool checkmesh()const;
private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif
	Real *points_;
	UInt *edges_;
	UInt *triangles_;
	UInt *neighbors_;

	UInt *border_edges; //contiene lista id_edges al bordo
	UInt num_nodes_, num_edges_, num_triangles_;

};


//!  SURFACE MESH:
//!  This class gives an object-oriented reading interface to the mesh object passed from R
/*!
 * The template parameters specify the order of its elements.
*/


template <UInt ORDER>
class MeshHandler<ORDER,2,3> {
public:
	typedef int UInt;
	//! A constructor.
    
    MeshHandler(Real* points, UInt* triangles,UInt* edges, UInt num_nodes, UInt num_triangles,UInt num_edges):
			points_(points), triangles_(triangles), edges_(edges), num_nodes_(num_nodes), num_triangles_(num_triangles),num_edges_(num_edges) {};
	
	//! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from a .csv file, useful for
      * debugging purposes
    */
	
    MeshHandler(std::string &filename){

       if(filename.find(".csv") != std::string::npos){
       		importfromCSV(filename);
       }
    }
	

    void importfromCSV(std::string &filename);
	
	//! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
    */
    #ifdef R_VERSION_
	MeshHandler(SEXP Rmesh);
	#endif

	~MeshHandler(){};

	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_nodes() const {return num_nodes_;}

	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */
    UInt num_triangles() const {return num_triangles_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning a Triangle
    /*!
     * \param id an Id argument
      \return The triangle with order coerent to that of the mesh with the specified id
    */
    Triangle<ORDER * 3,2,3>  getTriangle(Id id) const;

    void printPoints(std::ostream & out);
    void printTriangles(std::ostream & out);
   

     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a simply research between all triangle of the mesh
     * \param point the point we want to locate
      \return The triangle that contains the point
    */
    Triangle<ORDER * 3,2,3> findLocationNaive(Point point) const;


private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif

	std::vector<Real> points_;
	std::vector<UInt> triangles_;
	std::vector<UInt> edges_;


	UInt num_nodes_, num_triangles_, num_edges_;

};

#include "mesh_imp.h"

#endif
