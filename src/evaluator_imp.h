#ifndef __EVALUATOR_IMP_HPP__
#define __EVALUATOR_IMP_HPP__

template <UInt ORDER>
void Evaluator<ORDER,2,2>::eval(Real* X, Real *Y, UInt length, const Real *coef, UInt order, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,2> current_element;
	// std::vector<Triangle<3*ORDER> > starting_triangles; Problem with alignment not solved
	// by http://eigen.tuxfamily.org/dox-devel/group__TopicUnalignedArrayAssert.html
	// starting_triangles.resize(1);
	Element<Nodes,2,2> starting_element;


	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	starting_element = mesh_.getElement(0);
	for (int i = 0; i<length; ++i)
	{
		current_point = Point(X[i],Y[i]);
		//current_triangle = mesh_.findLocationNaive(current_point);
		//std::cout<<"Finding point.. "<<i<<" from "<<current_triangle.getId()<<"\n";
		current_element = mesh_.findLocationWalking(current_point, starting_element);
		//current_triangle.print(cout);
		//std::cout<<"Walking...triangle: "<< current_triangle.getId()<<std::endl;
		if(current_element.getId() == Identifier::NVAL && redundancy == true)
		{
			//To avoid problems with non convex mesh
			//std::cout<<"Position Not Found Walking... \n";
			current_element = mesh_.findLocationNaive(current_point);
			//std::cout<<"Naively...triangle: "<< current_triangle.getId()<<std::endl;
		}
		if(current_element.getId() == Identifier::NVAL)
		{
			//std::cout<<"Position Not Found Naively... \n";
			isinside[i]=false;
		}
		else
		{
			isinside[i]=true;
			for (int j=0; j<Nodes; ++j)
			{
				coefficients[j] = coef[current_element[j].getId()];
			}
			result[i] = evaluate_point<Nodes,2,2>(current_element, current_point, coefficients);
			starting_element = current_element;
		}
	}
}
;


template <UInt ORDER>
void Evaluator<ORDER,2,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, UInt order, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,3> current_element;
	Element<Nodes,2,3> starting_element;


	Point current_point;

	Eigen::Matrix<Real,Nodes,1> coefficients;
	starting_element = mesh_.getElement(0);
	for (int i = 0; i<length; ++i)
	{
		current_point = Point(X[i],Y[i],Z[i]);
		current_element = mesh_.findLocationNaive(current_point);

		if(current_element.getId() == Identifier::NVAL)
		{
			isinside[i]=false;
		}
		else
		{
			isinside[i]=true;
			for (int j=0; j<(Nodes); ++j)
			{
				coefficients[j] = coef[current_element[j].getId()];
			}
			result[i] = evaluate_point<Nodes,2,3>(current_element, current_point, coefficients);
			//std::cout<<"result = " <<result[i]<<"\n";
			starting_element = current_element;
		}
	}
}


template <UInt ORDER>
void Evaluator<ORDER,3,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, UInt order, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 6*ORDER-2;
	Element<Nodes,3,3> current_element;
	Element<Nodes,3,3> starting_element;


	Point current_point;

	Eigen::Matrix<Real,Nodes,1> coefficients;
	starting_element = mesh_.getElement(0);
	for (int i = 0; i<length; ++i)
	{
		current_point = Point(X[i],Y[i],Z[i]);
		current_element = mesh_.findLocationNaive(current_point);

		if(current_element.getId() == Identifier::NVAL)
		{
			isinside[i]=false;
		}
		else
		{
			isinside[i]=true;
			for (int j=0; j<Nodes; ++j)
			{
				coefficients[j] = coef[current_element[j].getId()];
			}
			result[i] = evaluate_point<Nodes,3,3>(current_element, current_point, coefficients);
			//std::cout<<"result = " <<result[i]<<"\n";
			starting_element = current_element;
		}
	}
}


#endif
