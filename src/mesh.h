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
			points_(points), edges_(edges), elements_(triangles), neighbors_(neighbors), num_nodes_(num_nodes), num_edges_(num_edges), num_elements_(num_triangles) {};

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
    UInt num_elements() const {return num_elements_;}

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

   //! A normal member setting an Element
        /*!
         * \param id an Id argument
          \return The element with order coerent to that of the mesh with the specified id
        */
    //void setElement(Element<3*ORDER,2,2>& tri, Id id) const;

    //! A normal member returning an Element
    /*!
     * \param id an Id argument 
      \return The element with order coerent to that of the mesh with the specified id
    */ 
    Element<3*ORDER,2,2>  getElement(Id id) const;

    //The "number" neighbor of element i is opposite the "number" corner of element i
    //! A normal member returning the Neighbors of a element
    /*!
     * \param id the id of the element
     * \param number the number of the vertex
      \return The element that has as an edge the one opposite to the specified
      vertex
    */
    Element<3*ORDER,2,2> getNeighbors(Id id_element, UInt number) const;

    void printPoints(std::ostream & out);
    void printEdges(std::ostream & out);
    void printElements(std::ostream & out);
    void printNeighbors(std::ostream & out);

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<3*ORDER,2,2> findLocationNaive(Point point) const;

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
     * \param point the point we want to locate
     * \param starting_elements a vector of points that specifies the poposed starting
     * points for the walking algorithm
      \return The element that contains the point
    */
    Element<3*ORDER,2,2> findLocationWalking(const Point& point, const Element<3*ORDER,2,2>& starting_element) const;

    //int readMesh(std::string const & file);
	//double measure()const;
	//bool checkmesh()const;
private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif
	Real *points_;
	UInt *edges_;
	UInt *elements_;
	UInt *neighbors_;

	UInt *border_edges; //contiene lista id_edges al bordo
	UInt num_nodes_, num_edges_, num_elements_;

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
    
    MeshHandler(Real* points, UInt* triangles, UInt num_nodes, UInt num_triangles):
			points_(points), elements_(triangles), num_nodes_(num_nodes), num_elements_(num_triangles) {};
	
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
    UInt num_elements() const {return num_elements_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Element
    /*!
     * \param id an Id argument
      \return The element with order coerent to that of the mesh with the specified id
    */
    Element<3*ORDER,2,3>  getElement(Id id) const;

    void printPoints(std::ostream & out);
    void printElements(std::ostream & out);
   

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<3*ORDER,2,3> findLocationNaive(Point point) const;


private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif

	std::vector<Real> points_;
	std::vector<UInt> elements_;


	UInt num_nodes_, num_elements_;

};


//!  VOLUME MESH:
//!  This class gives an object-oriented reading interface to the mesh object passed from R
/*!
 * The template parameters specify the order of its elements.
*/


template <UInt ORDER>
class MeshHandler<ORDER,3,3> {
public:
	typedef int UInt;
	//! A constructor.
    
    MeshHandler(Real* points, UInt* tetrahedrons, UInt num_nodes, UInt num_tetrahedrons):
			points_(points), elements_(tetrahedrons), num_nodes_(num_nodes), num_elements_(num_tetrahedrons) {};
	
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
      \return The number of elements in the mesh
    */
    UInt num_elements() const {return num_elements_;}

    //! A normal member returning a Point
    /*!
     * \param id an Id argument
      \return The point with the specified id
    */
    Point getPoint(Id id);

    //! A normal member returning an Element
    /*!
     * \param id an Id argument
      \return The element with order coerent to that of the mesh with the specified id
    */
    Element<6*ORDER-2,3,3>  getElement(Id id) const;

    void printPoints(std::ostream & out);
    void printElements(std::ostream & out);
   

     //! A normal member returning the element on which a point is located
    /*!
     * This method implements a simply research between all the elements of the mesh
     * \param point the point we want to locate
      \return The element that contains the point
    */
    Element<6*ORDER-2,3,3> findLocationNaive(Point point) const;


private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif

	std::vector<Real> points_;
	std::vector<UInt> elements_;


	UInt num_nodes_, num_elements_;

};





#include "mesh_imp.h"

#endif
