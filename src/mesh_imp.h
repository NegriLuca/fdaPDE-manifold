#ifndef MESH_IMP_H_
#define MESH_IMP_H_

#include<iostream>
#include<fstream>
#include<sstream>

#ifdef R_VERSION_
template <UInt ORDER>
MeshHandler<ORDER,2,2>::MeshHandler(SEXP mesh)
{
	mesh_ 		= mesh;
	points_ 	= REAL(VECTOR_ELT(mesh_, 0));
	edges_ 		= INTEGER(VECTOR_ELT(mesh_, 6));
	triangles_  = INTEGER(VECTOR_ELT(mesh_, 3));
	neighbors_  = INTEGER(VECTOR_ELT(mesh_, 8));

	num_nodes_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 0), R_DimSymbol))[0];
	num_edges_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 6), R_DimSymbol))[0];
	num_triangles_ = INTEGER(Rf_getAttrib(VECTOR_ELT(mesh_, 3), R_DimSymbol))[0];

}
#endif

template <UInt ORDER>
Point MeshHandler<ORDER,2,2>::getPoint(Id id)
{
	Point point(id, Identifier::NVAL, points_[id], points_[num_nodes_+id]);
	return point;
}

template <UInt ORDER>
Edge MeshHandler<ORDER,2,2>::getEdge(Id id)
{
	Id id_start_point = edges_[id];
	Id id_end_point = edges_[num_edges_+id];
	Edge edge(id, Identifier::NVAL, Point(id_start_point, Identifier::NVAL, points_[id_start_point], points_[num_nodes_+id_start_point]),
						Point(id_end_point, Identifier::NVAL, points_[id_end_point], points_[num_nodes_+id_end_point]));
	return edge;
}

template <UInt ORDER>
Triangle<ORDER * 3,2,2> MeshHandler<ORDER,2,2>::getTriangle(Id id) const
{
	std::vector<Point> triangle_points;
	triangle_points.resize(ORDER * 3);
	Id id_current_point;
	for (int i=0; i<ORDER * 3; ++i)
	{
		id_current_point = triangles_[i*num_triangles_ + id];
		triangle_points[i]= Point(id_current_point, Identifier::NVAL, points_[id_current_point],points_[num_nodes_+id_current_point]);
	}
	return Triangle<ORDER * 3,2,2>(id, triangle_points);
}

template <UInt ORDER>
Triangle<ORDER * 3,2,2>MeshHandler<ORDER,2,2>::getNeighbors(Id id_triangle, UInt number) const
{
	Id id_neighbour = neighbors_[number * num_triangles_ + id_triangle];
	//std::cout<<"Neighbour id "<< id_neighbour;
	if (id_neighbour == -1) return Triangle<ORDER * 3,2,2>(); //Triangle with NVAL ID

	return getTriangle(id_neighbour);
}

template <UInt ORDER>
Triangle<ORDER * 3,2,2> MeshHandler<ORDER,2,2>::findLocationNaive(Point point) const
{
	Triangle<ORDER * 3,2,2> current_triangle;
	//std::cout<<"Start searching point naively \n";
	for(Id id=0; id < num_triangles_; ++id)
	{
		current_triangle = getTriangle(id);
		if(current_triangle.isPointInside(point))
			return current_triangle;
	}
	//std::cout<<"Point not found \n";
	return Triangle<ORDER * 3,2,2>(); //default triangle with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER>
Triangle<ORDER * 3,2,2> MeshHandler<ORDER,2,2>::findLocationWalking(const Point& point, const Triangle<ORDER * 3,2,2>& starting_triangle) const
{

	//Walking algorithm to the point
	Triangle<ORDER * 3,2,2> current_triangle = starting_triangle;

	int direction=0;

	//Test for found Triangle, or out of border
	while(current_triangle.getId() != Identifier::NVAL && !current_triangle.isPointInside(point) )
	{
		direction = current_triangle.getPointDirection(point);
		//std::cout<<"Direction "<<direction<<";";
		current_triangle = getNeighbors(current_triangle.getId(), direction);
  	    //std::cout<<" ID "<<current_triangle.getId();
	}

	return current_triangle;
}


/*std::ostream & operator<<(std::ostream & out, MeshHandler const& m){
	out<< " ***** MESH  INFORMATION ******"<<std::endl;
	out<<" Num Points="<<m.num_nodes()<<" "<<" Num elements="<<m.num_triangles()<<" "
			<<"Num. edges="<<m.num_edges();
			//<<" "<<"Num Boundary Edges="<<m.num_bEdges()<<std::endl;
	out<< "POINTS:"<<std::endl;
	int oprec=out.precision(10);
	std::ios_base::fmtflags oflags=
			out.setf(std::ios_base::scientific,std::ios_base::floatfield);
	for (UInt i=0;i<m.num_nodes();++i){
		Point p=m.point(i);
		double x=p[0];
		double y=p[1];
		out<<i<<" "<<x<<" "<<y<<std::endl;
	}
	out<<" TRIANGLE CONNECTIVITY AND AREA:"<<std::endl;
	for (UInt i=0; i<m.num_elements();++i){
		Triangle t=m.triangle(i);
		out<<i<<" "<<t[0].id()<<" "<<t[1].id()<<" "<<t[2].id()<<
		  " "<<t.measure()<<"  Edge: "<<t.getEdges_id(0)<<" "<<t.getEdges_id(1)<<" "<<t.getEdges_id(2)<<std::endl;
	}
	out.precision(oprec);
	out.flags(oflags);
	return out;
}*/

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printPoints(std::ostream & out)
{
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[i]<<","<<points_[num_nodes_+i]<<")"<<std::endl<<"------"<<std::endl;
	}
}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printEdges(std::ostream & out)
{

	out << "Numero lati: "<< num_edges_ <<std::endl;
	for (UInt i = 0; i < num_edges_; ++i )
	{
		out<<"Lato ("<<edges_[i]<<","<<edges_[num_edges_+i]<<")"<<std::endl;
	}

}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printTriangles(std::ostream & out)
{

	out << "# Triangles: "<< num_triangles_ <<std::endl;
	for (UInt i = 0; i < num_triangles_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < ORDER * 3; ++k)
			out<<triangles_[k*num_triangles_ + i]<<"   ";
		out<<std::endl;
	}

}

template <UInt ORDER>
void MeshHandler<ORDER,2,2>::printNeighbors(std::ostream & out)
{

	out << "# Neighbors list: "<< num_triangles_ <<std::endl;
	for (UInt i = 0; i < num_triangles_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < 3; ++k)
			out<<neighbors_[k*num_triangles_ + i]<<"   ";
		out<<std::endl;
	}

}

//////////////////////////////////////////////////////////
// Implementation of class MeshHandler for surface mesh //
//////////////////////////////////////////////////////////


#ifdef R_VERSION_
template <UInt ORDER>
MeshHandler<ORDER,2,3>::MeshHandler(SEXP mesh)
{
	mesh_ = mesh;
	num_nodes_ = INTEGER(VECTOR_ELT(mesh_,0))[0];
	num_triangles_ = INTEGER(VECTOR_ELT(mesh_,1))[0];
	num_edges_ = INTEGER(VECTOR_ELT(mesh_,2))[0];
	points_.assign(REAL(VECTOR_ELT(mesh_, 3)) , REAL(VECTOR_ELT(mesh_, 3)) + 3*num_nodes_);
	triangles_.assign(INTEGER(VECTOR_ELT(mesh_, 4)), INTEGER(VECTOR_ELT(mesh_, 4))+ 3*ORDER*num_triangles_);
	edges_.assign(INTEGER(VECTOR_ELT(mesh_, 5)), INTEGER(VECTOR_ELT(mesh_, 5))+ 2*num_edges_);
	std::for_each(triangles_.begin(), triangles_.end(), [](int& i){i-=1;});
	std::for_each(edges_.begin(), edges_.end(), [](int& i){i-=1;});
}
#endif


template <UInt ORDER>
void MeshHandler<ORDER,2,3>::importfromCSV(std::string &filename){

	UInt nnodes;
	UInt ntriangles;
	UInt point_index;
	std::string line;
	std::string dummy;
	char comma;

	std::ifstream file;
	file.open(filename);


	// Read the number of points
	getline(file,line);
	std::istringstream ss(line);

	ss >> dummy; // throw away "num_points"
	ss >> nnodes;

	num_nodes_ = nnodes;
	points_.resize(3*nnodes);

	// Read the number of points
	getline(file,line);
	std::istringstream ss2(line);

	ss2 >> dummy; // throw away "num_triangles"
	ss2 >> ntriangles;

	num_triangles_ = ntriangles;
	triangles_.resize(3*ORDER*ntriangles);


	getline(file,line); //skip a white line

	// READ THE VERTICES MATRIX
	for(UInt i=0; i<nnodes; ++i){
		std::getline(file,line);
		std::istringstream ss(line);
		ss>>points_[3*i];
		ss>>comma;
		ss>>points_[3*i+1];
		ss>>comma;
		ss>>points_[3*i+2];
	};

	getline(file,line); //skip a white line

	// READ THE CONNECTIVIY MATRIX


	for(UInt i=0; i<ntriangles; ++i){
		std::getline(file,line);
		std::istringstream ss(line);
		ss>>point_index;
		triangles_[i*3] = --point_index;
		ss>>comma;
		ss>>point_index;
		triangles_[i*3+1] = --point_index;
		ss>>comma;
		ss>>point_index;
		triangles_[i*3+2] = --point_index;
		/*for(UInt k=0; k< 3*ORDER; ++k){
			ss>>point_index;
			triangles_[i*3*ORDER + k] = --point_index;
			ss>>comma;

		};*/

	};


};


template <UInt ORDER>
Point MeshHandler<ORDER,2,3>::getPoint(Id id)
{
	Point point(id, Identifier::NVAL, points_[id], points_[id+1],points_[id+2]);
	return point;
}

template <UInt ORDER>
Triangle<ORDER * 3,2,3> MeshHandler<ORDER,2,3>::getTriangle(Id id) const
{
	std::vector<Point> triangle_points;
	triangle_points.resize(ORDER * 3);
	Id id_current_point;
	for (int i=0; i<ORDER * 3; ++i)
	{
		id_current_point = triangles_[3*ORDER * id + i];
		triangle_points[i]= Point(id_current_point, Identifier::NVAL, points_[3*id_current_point],points_[3*id_current_point+1],points_[3*id_current_point+2]);
	}
	return Triangle<ORDER * 3,2,3>(id, triangle_points);
}

template <UInt ORDER>
Triangle<ORDER * 3,2,3> MeshHandler<ORDER,2,3>::findLocationNaive(Point point) const
{
	Triangle<ORDER * 3,2,3> current_triangle;
	//std::cout<<"Start searching point naively \n";
	for(Id id=0; id < num_triangles_; ++id)
	{
		current_triangle = getTriangle(id);
		if(current_triangle.isPointInside(point))
			return current_triangle;
	}
	//std::cout<<"Point not found \n";
	return Triangle<ORDER * 3,2,3>(); //default triangle with NVAL ID
}

template <UInt ORDER>
void MeshHandler<ORDER,2,3>::printPoints(std::ostream & out)
{
std::cout<<"printing points"<<"\n";
	for(UInt i = 0; i < num_nodes_; ++i)
	{
		out<<"-"<< i <<"-"<<"("<<points_[3*i]<<","<<points_[3*i+1]<<","<<points_[3*i+2]<<")"<<std::endl<<"------"<<std::endl;
	}
}

template <UInt ORDER>
void MeshHandler<ORDER,2,3>::printTriangles(std::ostream & out)
{

	out << "# Triangles: "<< num_triangles_ <<std::endl;
	for (UInt i = 0; i < num_triangles_; ++i )
	{
		out<<"-"<< i <<"- ";
		for( UInt k = 0; k < ORDER * 3; ++k)
			out<<triangles_[i*3*ORDER + k]<<"   ";
		out<<std::endl;
	}

}


#endif
