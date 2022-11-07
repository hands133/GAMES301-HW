#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

#include <tuple>
#include <vector>
#include <functional>
#include <Eigen\Dense>
#include <Eigen\Sparse>

#include <TinyAD\Scalar.hh>

// 2D, Symmetric Dirichlet Energy Only
namespace eigensys
{
	struct QPW_Invariables
	{
		void CalculateInvariants(const Eigen::Matrix2d& S);
		// Initialized with unit matrix S
		double I1{ 2.0 };	// I1 = sum(sigma_i)
		double I2{ 2.0 };	// I2 = sum(sigma_i^2)
		double I3{ 1.0 };	// I3 = prod(sigma_i)
	};

	struct QPW_Decomposition
	{
		void CalculateSVDnPolar(const Eigen::Matrix2d& F);
		Eigen::Matrix2d U, Sigma, V;	// { U, Sigma, V } = SVD(F)
		Eigen::Matrix2d S, R;			// polar decomposition
	};

	struct QPW_DeformGradient
	{
		void CalculateDeformGradient(const Eigen::Matrix2d& Ds, const Eigen::Matrix2d& DmInv);
		Eigen::Matrix2d F{ Eigen::Matrix2d::Identity() };	// deformation gradient F
	};

	//struct QPW_EnergySD_2DGradient;

	struct QPW_DeformVectors
	{
		void CalculateVectors(const Eigen::Matrix2d& F, const QPW_Decomposition& decomp);
		
		Eigen::Vector4d f{ Eigen::Vector4d::Zero() };		// vec(F)
		Eigen::Vector4d g{ Eigen::Vector4d::Zero() };		// vec(R)
		Eigen::Vector4d r{ Eigen::Vector4d::Zero() };		// r / sqrt(2)
		Eigen::Vector4d t{ Eigen::Vector4d::Zero() };		// vec(U twist VT) / sqrt(2)
		Eigen::Vector4d p{ Eigen::Vector4d::Zero() };		// vec(U [1 -1] VT) / sqrt(2)
		Eigen::Vector4d l{ Eigen::Vector4d::Zero() };		// vec(U {1  1} VT) / sqrt(2)
		Eigen::Vector4d d1{ Eigen::Vector4d::Zero() };		// vec(U {1  1} VT) / sqrt(2)
		Eigen::Vector4d d2{ Eigen::Vector4d::Zero() };		// vec(U {1  1} VT) / sqrt(2)
	};
	
	// F Dm = Ds
	class QPW_DataPack
	{
		//friend struct QPW_pPSIq_pfq;
		friend struct QPW_EnergySD_2D;
		//friend struct QPW_EnergySD_2DGradient;
		friend class ProjectNewtonSolver;
	public:
		QPW_DataPack(const Eigen::Matrix2d& Dm);

		void CalculatepF_pDs(const Eigen::Matrix2d& Ds);

		// call CalculateGradient ahead!
		Eigen::Matrix<double, 4, 6> GetLocalpfq_pxq() const { return pfq_pxq; }
		Eigen::MatrixXd GetGlobalpf_px(acamcad::polymesh::PolyMesh* mesh, size_t faceID) const;

	private:
		void Calculatepfq_pxq(const Eigen::Matrix2d& DmInv);

	private:
		Eigen::Matrix2d DmINV;	// F = Ds / Dm
		Eigen::Matrix<double, 4, 6> pfq_pxq{ Eigen::Matrix<double, 4, 6>::Zero() };	// pf / px per quad-point

		QPW_DeformGradient	m_DeformGradient;
		QPW_Decomposition	m_Decomposition;
		QPW_Invariables		m_Invariables;
		QPW_DeformVectors	m_DeformVectors;
	};

	//struct QPW_pPSIq_pfq
	//{
	//	Eigen::Vector4d operator()(const QPW_pfq_pxq& pfq_pxq);
	//};

	//// quadrature-point-wise Energy Gradient
	//struct QPW_EnergySD_2DGradient
	//{
	//	void CalculateEnergyGradient(const QPW_pfq_pxq& pfq_pxq);
	//	Eigen::Vector<double, 6> gradient = Eigen::Vector<double, 6>::Zero();
	//	Eigen::Matrix<double, 6, 6> hessian = Eigen::Matrix<double, 6, 6>::Zero();

	//private:
	//	Eigen::Vector<double, 6> CalculateGrad(const QPW_pfq_pxq& pfq_pxq) const;
	//	Eigen::Matrix<double, 6, 6> CalculateHess(const QPW_pfq_pxq& pfq_pxq) const;
	//};

	// global Energy Gradient
	struct EnergySD_2DGradient
	{
		Eigen::MatrixXd operator()(
			acamcad::polymesh::PolyMesh* mesh,
			const Eigen::VectorXd& UVs, 
			std::vector<QPW_DataPack>& packs);
	};

	struct QPW_EigenValVecPair
	{
		double l{ 0.0 };
		Eigen::Vector4d e{ Eigen::Vector4d::Zero() };
	};

	struct QPW_EigenSystem2D
	{
		std::array<QPW_EigenValVecPair, 4> valvecpairs;
	};

	class ProjectNewtonSolver
	{
	public:
		void PresetMeshUV(acamcad::polymesh::PolyMesh* mesh);
		// Update Mesh Once a time
		bool UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh);

	private:
		Eigen::VectorXd	SetupUVs(acamcad::polymesh::PolyMesh* mesh);
		
		QPW_EigenSystem2D Eval_Energy_EigenSystem(const QPW_DataPack& pfq_pxq) const;
		
		Eigen::Matrix4d Calculatep2PSIq_pfq2(const QPW_EigenSystem2D& eigensys) const;
		Eigen::Vector4d CalculatepPSIq_pfq(const QPW_DataPack& pack);

		std::pair<double, double> Line_Search(
			acamcad::polymesh::PolyMesh* polymesh,
			const Eigen::VectorXd& d,
			const Eigen::VectorXd& grad,
			const Eigen::VectorXd& UVs,
			double lastEnergy,
			double gamma, double c);

		double CalculateQPWEnergySD_2DNoArea(const Eigen::Matrix2d& DmINV, const Eigen::Matrix2d& Ds) const;
		double CalculateQPWEnergySD_2DNoArea(const QPW_Invariables& invars) const;

		double CalculateEnergySD_2D(
			acamcad::polymesh::PolyMesh* polymesh,
			const Eigen::VectorXd& UVs,
			const std::vector<Eigen::Matrix2d>& DsList) const;

		std::tuple<Eigen::VectorXd, Eigen::SparseMatrix<double>> CalculateGlobalEnergyDerivative(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs);
		std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> CalculateLocalEnergyDerivativeNoArea(
			acamcad::polymesh::PolyMesh* mesh, 
			const Eigen::Matrix2d& Dm,
			const Eigen::Matrix2d& Ds,
			const Eigen::VectorXd& UVs,
			size_t faceID);

	private:
		//std::vector<QPW_pfq_pxq> m_pf_pxList;
		std::vector<Eigen::Matrix2d> m_DmList;
		Eigen::VectorXd m_UVList;

		size_t Iters = 0;
		double m_Energy = 0.0;

		bool m_FirstUpdate = true;
	};
}