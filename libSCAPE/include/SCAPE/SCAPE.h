#pragma once

#include "MatrixnD.h"
#include "SCAPE_structs.h"
#include "PCA_basis.h"
#include <umfpack/umfpack.h>
#include "operators/alignPointClouds.h"
#include "RandAccSparseMatrix.h"
#include <omp.h>

//#define SAVE_MATRIX 0

typedef SquareMatrixND<vec3d> SMat3D;

static const int numParts = 16;
static const int numJoints = 15;
static const int numT = 25000;

//! Can be used to store the rigid per part rotations.
struct PartRotations {
	SMat3D R[numParts];
};

//! Can be used to store the rigid per joint rotation.
struct JointRotations {
	vec3d rv[numJoints];
};

class SCAPE_Model {

public:
	SCAPE_Model() {
		NumericQ = NULL;
		SymbolicQ = NULL;
		NumericDGrad = NULL;
		SymbolicDGrad = NULL;
	}

	~SCAPE_Model() {
		if (NumericQ) umfpack_di_free_numeric(&NumericQ);
		if (SymbolicQ) umfpack_di_free_symbolic(&SymbolicQ);
		if (NumericDGrad) umfpack_di_free_numeric(&NumericDGrad);
		if (SymbolicDGrad) umfpack_di_free_symbolic(&SymbolicDGrad);
	}


	//! Returns the two joints that define the non-rigid pose deformation of the triangle 'id'. If there is only one joint, the second joint ID will be set to -1!
	vec2i getClosestJointsForTriangle(UInt id) const {
		return triToJoints[id];
	}


	//! Returns the non-rigid pose deformation of the triangle 'id' for joint rotations 'j0_rot' and 'j1_rot'.
	SMat3D estimateQfromJoints(UInt id, const vec3d& j0_rot, const vec3d& j1_rot) const {
		SMat3D Qret;
		for (int x = 0; x < 3; x++) {
			for (int y = 0; y < 3; y++) {
				const double* a = perTriangleParameters[id].vals[x][y];
				Qret(x, y) = j0_rot.x * a[0] + j0_rot.y * a[1] + j0_rot.z * a[2] + j1_rot.x * a[3] + j1_rot.y * a[4] + j1_rot.z * a[5] + a[6];
			}
		}
		return Qret;
	}

	RigidTransform<vec3d> computePartAlignment(const Part& part, const std::vector<vec3d>& template_verts, const std::vector<vec3d>& pose_verts) const {

		std::vector<vec3d> pts_template;
		std::vector<vec3d> pts_pose;
		for (UInt i1 = 0; i1 < part.verts.size(); i1++) {
			pts_template.push_back(template_verts[part.verts[i1]]);
			pts_pose.push_back(pose_verts[part.verts[i1]]);
		}

		RigidTransform<vec3d> trans = rigidAlignPointClouds(&pts_pose[0], &pts_template[0], pts_template.size());
		return trans;
	}


	std::vector<double> projectDToPCABasis(std::vector<SMat3D>& D, int cropToNBasisVectors = 0) const {
		std::vector<double> dlin(numT * 9);
		double* Dpoint = (double*)(&D[0]);
		for (int i1 = 0; i1 < numT * 9; i1++) {
			dlin[i1] = Dpoint[i1] - avg_defgrad[i1];
		}
		std::vector<double> inpca = pcabasis.projectIntoBasis(dlin);

		if (cropToNBasisVectors != 0) {
			for (UInt i1 = cropToNBasisVectors; i1 < inpca.size(); i1++) inpca[i1] = 0;
		}

		std::vector<double> reproj = pcabasis.projectFromBasis(inpca);
		for (int i1 = 0; i1 < numT * 9; i1++) {
			reproj[i1] += avg_defgrad[i1];
		}

		for (int i1 = 0; i1 < numT * 9; i1++) {
			Dpoint[i1] = reproj[i1];
		}

		return inpca;
	}


	/// <summary>
	/// Create a new instance with a given pose and shape.
	/// 
	/// The pose consists of part rotations R_l and per-triangle deformations Q_k.
	/// The shape consists of per-triangle deformations S_k.
	/// 
	/// The preliminary resulting vertices are given by equation (6) in [1]:
	/// 
	///   v_{k, j} = R_l[k] S_k Q_k \hat{v}_{k, j}
	/// 
	/// where k ranges over the set of triangles and j over {1, 2, 3} (three points per triangle).
	/// 
	/// TODO: Really R_l[k] and not part rotations?
	/// 
	/// Note that for a complete consistent mesh multiple triangles must share common vertices,
	/// i.e. they must agree upon them. Therefore, the final vertices are determined in such a way
	/// that they minimize the compromise needed by every triangle, cf. equation (9) in [1].
	/// </summary>
	/// <param name="pr">The part rotations R_l for each part l.</param>
	/// <param name="Q">The 3x3 linear transformation matrices Q_k deforming every triangle k.</param>
	/// <param name="D">The 3x3 shape deformation matrices S_k deforming every triangle k.</param>
	/// <returns>TODO</returns>
	/// 
	/// <seealso cref="decodeDefGradMesh"/>
	std::vector<vec3d> reconstructModel(const PartRotations* pr, const std::vector<SMat3D>* Q, const std::vector<SMat3D>* D) const {
		UInt num;
		if (Q != NULL) num = Q->size();
		else num = D->size();

		std::vector< SquareMatrixND<vec3d> > Rots(num);

		if (pr == NULL || Q == NULL) { // Only body shape change
			std::cerr << "Only shape" << std::endl;
#pragma omp parallel
			for (UInt i2 = 0; i2 < numT; i2++) Rots[i2] = (*D)[i2];
		} else if (D == NULL) { // default body shape in input pose
#pragma omp parallel
			// TODO: Why part rotations and not join rotations as per equation (6) in [1]
			for (UInt i2 = 0; i2 < numT; i2++) Rots[i2] = pr->R[triIDtoPartID[i2]] * (*Q)[i2];
		} else { // transfer pose and shape
#pragma omp parallel
			for (UInt i2 = 0; i2 < numT; i2++) Rots[i2] = (*D)[i2] * pr->R[triIDtoPartID[i2]] * (*Q)[i2];
		}


		return decodeDefGradMesh(Rots);
	}

	//! Computes and returns the rigid rotations of the mesh parts.
	PartRotations getRigidPartRotations(const std::vector<vec3d>& pose_verts) const {

		PartRotations pr;
		for (UInt i1 = 0; i1 < numParts; i1++) {
			RigidTransform<vec3d> trans = computePartAlignment(rigidParts[i1], default_mesh_vertices, pose_verts);
			pr.R[i1] = trans.rotate;
		}
		return pr;

	}


	//! Returns the joint rotations as rotation vectors. Note that the rotation vectors are already projected into their PCA basis.
	JointRotations getJointRotationsFromPartRotations(const PartRotations& pr) const {
		JointRotations jr;
		// Compute joint rotations
		for (UInt i1 = 0; i1 < numJoints; i1++) {
			const Joint& j = joints[i1];
			SquareMatrixND<vec3d> RM = pr.R[j.part1].getTransposed() * pr.R[j.part0];
			jr.rv[i1] = RM.getRotVec();
			jr.rv[i1] *= -1;

			// Proj to pca basis
			jr.rv[i1] = joint_PCA_bases[i1].TransposedVecTrans(jr.rv[i1]);
		}
		return jr;
	}

	std::vector<SMat3D> getAllQs(const JointRotations& jr) const {

		std::vector<SMat3D> Qs(default_mesh_triangles.size());

#pragma omp parallel
		for (UInt i1 = 0; i1 < default_mesh_triangles.size(); i1++) {
			vec2i js = triToJoints[i1];
			const vec3d& j0_rot = jr.rv[js.x];
			vec3d j1_rot(0, 0, 0);
			if (js.y > -1) j1_rot = jr.rv[js.y];
			SMat3D Q = estimateQfromJoints(i1, j0_rot, j1_rot);
			Qs[i1] = Q;
		}

		return Qs;
	}

	//! Projects the input model to SCAPE Space.
	std::vector<vec3d> projectToSCAPESpace(const std::vector<vec3d>& verts, bool poseOnly = false, const std::vector<double>& shapeParams = std::vector<double>()) const {
		PartRotations pr = getRigidPartRotations(verts);
		JointRotations jr = getJointRotationsFromPartRotations(pr);

		// pose
		std::vector<SquareMatrixND<vec3d>> Q = getAllQs(jr);
		std::vector<SquareMatrixND<vec3d>> Rots;
		for (UInt i2 = 0; i2 < numT; i2++) {
			SquareMatrixND<vec3d> R = pr.R[triIDtoPartID[i2]] * Q[i2];
			Rots.push_back(R);
		}

		// shape
		std::vector<SquareMatrixND<vec3d>> D;
		computeDMatrices(verts, Rots, D);
		std::vector<double> inpca = projectDToPCABasis(D, 0);

		// override shape params
		if (shapeParams.size() == pcabasis.numVecs) {
			double* Dp = (double*)&D[0];
			std::vector<double> Dlin = pcabasis.projectFromBasis(shapeParams);
#pragma omp parallel
			for (int i1 = 0; i1 < numT * 9; i1++) Dp[i1] = Dlin[i1] + avg_defgrad[i1];
		}

		if (poseOnly) {
			return reconstructModel(&pr, &Q, NULL); // Reconstruct pose only
		} else {
			return reconstructModel(&pr, &Q, &D); // Reconstruct pose and shape
		}
	}

	//! Projects the input model to SCAPE Space.
	std::vector<vec3d> projectToSCAPESpace(const std::vector<vec3d>& verts, std::vector<double>* params = NULL, int crop = 0) const {
		PartRotations pr = getRigidPartRotations(verts);
		JointRotations jr = getJointRotationsFromPartRotations(pr);

		std::vector<SquareMatrixND<vec3d>> Q = getAllQs(jr);
		std::vector<SquareMatrixND<vec3d>> Rots;
		for (UInt i2 = 0; i2 < numT; i2++) {
			SquareMatrixND<vec3d> R = pr.R[triIDtoPartID[i2]] * Q[i2];
			Rots.push_back(R);
		}
		std::vector<SquareMatrixND<vec3d>> D;
		computeDMatrices(verts, Rots, D);

		std::vector<double> inpca = projectDToPCABasis(D, crop);

		if (params != NULL) {
			params->resize(3 * numParts + pcabasis.numVecs);
			for (int i1 = 0; i1 < numParts; i1++) {
				vec3d partRot = pr.R[i1].getRotVec();
				(*params)[i1 * 3 + 0] = partRot.x;
				(*params)[i1 * 3 + 1] = partRot.y;
				(*params)[i1 * 3 + 2] = partRot.z;
			}
			for (int i1 = 0; i1 < pcabasis.numVecs; i1++) {
				(*params)[3 * numParts + i1] = inpca[i1];
			}
		}

		return reconstructModel(&pr, &Q, &D); // Reconstruct pose and shape
	}


	/// <summary>
	/// Create a vertices model from a (pose|shape) vector.
	/// </summary>
	/// <param name="params">
	/// A combined (pose|shape) vector where pose = (twist_1, twist_2, ..., twist_{numParts}) and
	/// shape = (pcaCoeff_1, ..., pcaCoeff_{pcabasis.numVecs}).
	/// 
	/// The pose subvector shall contain <see cref="numParts"/> many three-dimensional twist subvectors for every part rotation.
	/// A twist vector represents the rotation in an axis-angle fashion: The axis is determined by the vector's direction and
	/// the angle is determined by the vector's magnitude.
	/// 
	/// The shape subvector shall contain <see cref="pcabasis.numvecs"/> many scalar PCA coefficients.
	/// </param>
	/// <returns>
	/// Vertices in such an order that combining them with <see cref="topologically_correct_mesh_triangles"/>
	/// leads to the desired triangle mesh model.
	/// </returns>
	std::vector<vec3d> reconstructFromSCAPEParams(const std::vector<double>& params) const {

		PartRotations pr;
		for (int i1 = 0; i1 < numParts; i1++) {
			vec3d rv(params[i1 * 3 + 0], params[i1 * 3 + 1], params[i1 * 3 + 2]);
			pr.R[i1].setToRotationMatrixNew(rv);
		}
		JointRotations jr = getJointRotationsFromPartRotations(pr);


		std::vector<SquareMatrixND<vec3d>> Q = getAllQs(jr);


		std::vector<SquareMatrixND<vec3d>> Rots(Q.size());

#pragma omp parallel
		for (UInt i2 = 0; i2 < numT; i2++) {
			SquareMatrixND<vec3d> R = pr.R[triIDtoPartID[i2]] * Q[i2];
			Rots[i2] = R;
		}
		std::vector<SquareMatrixND<vec3d>> D(numT);
		double* Dp = (double*)&D[0];

		std::vector<double> inpca(pcabasis.numVecs);

		for (int i1 = 0; i1 < pcabasis.numVecs; i1++) inpca[i1] = params[3 * numParts + i1];
		std::vector<double> Dlin = pcabasis.projectFromBasis(inpca);

#pragma omp parallel
		for (int i1 = 0; i1 < numT * 9; i1++) Dp[i1] = Dlin[i1] + avg_defgrad[i1];

		return reconstructModel(&pr, &Q, &D); // Reconstruct pose and shape
	}


	//! Rotates the model to match the orientation of the base mesh
	RigidTransform<vec3d> rigidlyAlignModel(std::vector<vec3d>& verts) const {
		RigidTransform<vec3d> rtrans = rigidAlignPointClouds(&default_mesh_vertices[0], &verts[0], verts.size());
		for (UInt i1 = 0; i1 < verts.size(); i1++) verts[i1] = rtrans.vecTrans(verts[i1]);
		return rtrans;
	}


	//! Rotates the model to match the orientation of the base mesh
	void rigidlyAlignModel(std::vector<vec3d>& verts, const std::vector<vec3d>& verts_reference, std::vector<double>* weights = NULL) const {
		RigidTransform<vec3d> rtrans;
		if (weights == NULL) rtrans = rigidAlignPointClouds(&verts_reference[0], &verts[0], verts.size());
		else rtrans = rigidAlignPointCloudsWeighted(&verts_reference[0], &verts[0], verts.size(), &(*weights)[0]);
		for (UInt i1 = 0; i1 < verts.size(); i1++) verts[i1] = rtrans.vecTrans(verts[i1]);
	}


	void saveToFile(const char* fn) const {
		std::ofstream fout(fn, std::ios::binary);
		fout << "scape model container   ";

		int numV = (int)default_mesh_vertices.size();
		fout.write((char*)&numV, sizeof(int));
		fout.write((char*)&default_mesh_vertices[0], sizeof(vec3d) * numV);

		int numT = (int)default_mesh_triangles.size();
		fout.write((char*)&numT, sizeof(int));
		fout.write((char*)&default_mesh_triangles[0], sizeof(vec3i) * numT);

		fout.write((char*)&numParts, sizeof(int));
		for (int i1 = 0; i1 < numParts; i1++) rigidParts[i1].savePart(fout);

		if (triToJoints.size() != numT) {
			std::cerr << "size of triToJoints doesn't match number of triangles!" << std::endl;
			std::cerr << triToJoints.size() << std::endl;
			std::cerr << numT << std::endl;
			exit(1);
		}
		fout.write((char*)&numT, sizeof(int));
		fout.write((char*)&triToJoints[0], sizeof(vec2i) * numT);

		if (perTriangleParameters.size() != numT) {
			std::cerr << "size of perTriangleParams doesn't match number of triangles!" << std::endl;
			exit(1);
		}
		fout.write((char*)&numT, sizeof(int));
		fout.write((char*)&perTriangleParameters[0], sizeof(perTriangleParams) * numT);
		std::cerr << "sizeof(perTriangleParams): " << sizeof(perTriangleParams) << std::endl;

		if (joint_PCA_bases.size() != numJoints) {
			std::cerr << "size of joint_PCA_bases doesn't match number of joints!" << std::endl;
			exit(1);
		}
		fout.write((char*)&numJoints, sizeof(int));
		fout.write((char*)&joint_PCA_bases[0], sizeof(SMat3D) * numJoints);
		std::cerr << "sizeof(SMat3D): " << sizeof(SMat3D) << std::endl;

		if (joints.size() != numJoints) {
			std::cerr << "size of joints doesn't match number of joints!" << std::endl;
			exit(1);
		}
		fout.write((char*)&numJoints, sizeof(int));
		fout.write((char*)&joints[0], sizeof(Joint) * numJoints);
		std::cerr << "sizeof(Joint): " << sizeof(Joint) << std::endl;

		if (triIDtoPartID.size() != numT) {
			std::cerr << "size of triIDtoPartID doesn't match number of triangles!" << std::endl;
			exit(1);
		}
		fout.write((char*)&numT, sizeof(int));
		fout.write((char*)&triIDtoPartID[0], sizeof(int) * numT);

		int numX = (int)avg_defgrad.size();
		fout.write((char*)&numX, sizeof(int));
		fout.write((char*)&avg_defgrad[0], numX * sizeof(double));

		// Save tri_neighbors:
		numX = tri_neighbors.size();
		fout.write((char*)&numX, sizeof(int));
		fout.write((char*)&tri_neighbors[0], sizeof(vec3i) * numX);

		pcabasis.saveToBinaryStream(fout);

		int numTcorrect = (int)topologically_correct_mesh_triangles.size();
		fout.write((char*)&numTcorrect, sizeof(int));
		fout.write((char*)&topologically_correct_mesh_triangles[0], sizeof(vec3i) * numTcorrect);

		numX = (int)avg_jointangles.size();
		fout.write((char*)&numX, sizeof(int));
		fout.write((char*)&avg_jointangles[0], numX * sizeof(double));
		pose_pcabasis.saveToBinaryStream(fout);

		numX = (int)avg_SCAPE_coeffs.size();
		fout.write((char*)&numX, sizeof(int));
		fout.write((char*)&avg_SCAPE_coeffs[0], numX * sizeof(double));

		fout.close();
	}

	//
	// init code
	//
	void readFromFile(const char* fn, const char* pathToMatrices, bool outputStatusMessagesToStderr = false) {

		std::ifstream fin(fn, std::ios::binary);
		if (!fin.good()) {
			std::cerr << "Error reading file " << fn << " in SCAPE_Model.readFromFile()" << std::endl;
			exit(1);
		}

		char header[24];
		fin.read(header, sizeof(char) * 24);

		//std::cerr << "SCAPE Header: " << header << std::endl;

		int numV;
		fin.read((char*)&numV, sizeof(int));
		if (outputStatusMessagesToStderr) {
			std::cerr << "numV: " << numV << std::endl;
		}
		default_mesh_vertices.resize(numV);
		fin.read((char*)&default_mesh_vertices[0], sizeof(vec3d) * numV);

		int numT;
		fin.read((char*)&numT, sizeof(int));
		if (outputStatusMessagesToStderr) {
			std::cerr << "numT: " << numT << std::endl;
		}
		default_mesh_triangles.resize(numT);
		fin.read((char*)&default_mesh_triangles[0], sizeof(vec3i) * numT);

		int numP;
		fin.read((char*)&numP, sizeof(int));
		if (numP != numParts) {
			std::cerr << "number of parts doesn't match numParts" << std::endl;
			exit(1);
		}
		for (int i1 = 0; i1 < numParts; i1++) rigidParts[i1].readPart(fin);

		int numTTJ;
		fin.read((char*)&numTTJ, sizeof(int));
		if (numTTJ != numT) {
			std::cerr << "numTTJ doesn't match numTriangles" << std::endl;
			exit(1);
		}
		triToJoints.resize(numT);
		fin.read((char*)&triToJoints[0], sizeof(vec2i) * numT);


		int numTP;
		fin.read((char*)&numTP, sizeof(int));
		if (numTP != numT) {
			std::cerr << "numTP doesn't match numTriangles" << std::endl;
			exit(1);
		}
		perTriangleParameters.resize(numT);
		fin.read((char*)&perTriangleParameters[0], sizeof(perTriangleParams) * numT);

		int numJB;
		fin.read((char*)&numJB, sizeof(int));
		if (numJB != numJoints) {
			std::cerr << "numJB doesn't match numJoints" << std::endl;
			exit(1);
		}
		joint_PCA_bases.resize(numJoints);
		fin.read((char*)&joint_PCA_bases[0], sizeof(SMat3D) * numJoints);

		int numJ;
		fin.read((char*)&numJ, sizeof(int));
		if (numJ != numJoints) {
			std::cerr << "numJ doesn't match numJoints" << std::endl;
			exit(1);
		}
		joints.resize(numJoints);
		fin.read((char*)&joints[0], sizeof(Joint) * numJoints);

		int numTTtP;
		fin.read((char*)&numTTtP, sizeof(int));
		if (numTTtP != numT) {
			std::cerr << "numTTtP doesn't match numTriangles" << std::endl;
			exit(1);
		}
		triIDtoPartID.resize(numT);
		fin.read((char*)&triIDtoPartID[0], sizeof(int) * numT);

		int numX;
		fin.read((char*)&numX, sizeof(int));
		avg_defgrad.resize(numX);
		fin.read((char*)&avg_defgrad[0], numX * sizeof(double));

		// Read tri_neighbors:
		int numTriN;
		fin.read((char*)&numTriN, sizeof(int));
		tri_neighbors.resize(numTriN);
		fin.read((char*)&tri_neighbors[0], sizeof(vec3i) * numTriN);

		// Read PCA Basis
		int rnv = 0;
		pcabasis.loadFromStream(fin, rnv);

		// Read manifold topology
		int numTcorrect;
		fin.read((char*)&numTcorrect, sizeof(int));
		topologically_correct_mesh_triangles.resize(numTcorrect);
		fin.read((char*)&topologically_correct_mesh_triangles[0], sizeof(vec3i) * numTcorrect);

		// Read pose pca basis
		if (fin.peek() != EOF) {
			if (outputStatusMessagesToStderr) {
				std::cerr << "Reading pose PCA basis" << std::endl;
			}

			fin.read((char*)&numX, sizeof(int));
			avg_jointangles.resize(numX);
			fin.read((char*)&avg_jointangles[0], numX * sizeof(double));

			int rnv = 0;
			pose_pcabasis.loadFromStream(fin, rnv);
		}

		// Read pose pca basis
		if (fin.peek() != EOF) {
			if (outputStatusMessagesToStderr) {
				std::cerr << "Reading avg coeffs" << std::endl;
			}

			fin.read((char*)&numX, sizeof(int));
			avg_SCAPE_coeffs.resize(numX);
			fin.read((char*)&avg_SCAPE_coeffs[0], numX * sizeof(double));
		}

		fin.close();

		loadUMFPACKMatrices(pathToMatrices);
	}


	void computeDMatrices(const std::vector<vec3d>& pose_verts, std::vector< SquareMatrixND<vec3d> > Rots, std::vector< SquareMatrixND<vec3d> >& D) const {

		//double w = 0.000000001;
		double w = 0.00001;

		D.clear();
		D.resize(numT);

		std::vector<EdgePair> edges_template(numT);
		std::vector<EdgePair> edges_pose(numT);

		for (UInt i1 = 0; i1 < numT; i1++) {
			const vec3i& t = default_mesh_triangles[i1];
			vec3d v0_template = pose_verts[t.y] - pose_verts[t.x];
			vec3d v1_template = pose_verts[t.z] - pose_verts[t.x];
			vec3d v0_pose = default_mesh_vertices[t.y] - default_mesh_vertices[t.x];
			vec3d v1_pose = default_mesh_vertices[t.z] - default_mesh_vertices[t.x];
			const SquareMatrixND<vec3d>& M = Rots[i1];
			edges_template[i1] = EdgePair(M.vecTrans(v0_pose), M.vecTrans(v1_pose));
			edges_pose[i1] = EdgePair(v0_template, v1_template);
		}

		// Build right hand side
		RandAccessCompRowMatrix<double> A_123(3 * numT, 3 * numT);
		std::vector<double> rhs_123[3];
		rhs_123[0].resize(3 * numT);
		rhs_123[1].resize(3 * numT);
		rhs_123[2].resize(3 * numT);

		//#pragma omp parallel
		for (UInt i1 = 0; i1 < numT; i1++) {

			// Insert the three lines for a row of the i-th triangle
			const EdgePair& ep0 = edges_template[i1];
			const EdgePair& ep1 = edges_pose[i1];

			SquareMatrixND<vec3d> mat(ep0.v0);
			mat.addFromTensorProduct(ep0.v1, ep0.v1);

			for (int x = 0; x < 3; x++) {
				for (int y = 0; y < 3; y++) {
					A_123.add(3 * i1 + x, 3 * i1 + y, mat(x, y));
				}
			}

			rhs_123[0][3 * i1 + 0] = ep0.v0.x * ep1.v0.x + ep0.v1.x * ep1.v1.x;
			rhs_123[0][3 * i1 + 1] = ep0.v0.y * ep1.v0.x + ep0.v1.y * ep1.v1.x;
			rhs_123[0][3 * i1 + 2] = ep0.v0.z * ep1.v0.x + ep0.v1.z * ep1.v1.x;
			rhs_123[1][3 * i1 + 0] = ep0.v0.x * ep1.v0.y + ep0.v1.x * ep1.v1.y;
			rhs_123[1][3 * i1 + 1] = ep0.v0.y * ep1.v0.y + ep0.v1.y * ep1.v1.y;
			rhs_123[1][3 * i1 + 2] = ep0.v0.z * ep1.v0.y + ep0.v1.z * ep1.v1.y;
			rhs_123[2][3 * i1 + 0] = ep0.v0.x * ep1.v0.z + ep0.v1.x * ep1.v1.z;
			rhs_123[2][3 * i1 + 1] = ep0.v0.y * ep1.v0.z + ep0.v1.y * ep1.v1.z;
			rhs_123[2][3 * i1 + 2] = ep0.v0.z * ep1.v0.z + ep0.v1.z * ep1.v1.z;

			// Add constraints:
			for (UInt x = 0; x < 3; x++) {
				const int& other = tri_neighbors[i1][x];
				A_123.add(3 * i1 + 0, 3 * i1 + 0, w); A_123.add(3 * i1 + 0, 3 * other + 0, -w);
				A_123.add(3 * i1 + 1, 3 * i1 + 1, w); A_123.add(3 * i1 + 1, 3 * other + 1, -w);
				A_123.add(3 * i1 + 2, 3 * i1 + 2, w); A_123.add(3 * i1 + 2, 3 * other + 2, -w);
			}

		}


		std::vector<double> entries;
		std::vector<int> row_index;
		std::vector<int> col_ptr;
		A_123.getMatrix(entries, row_index, col_ptr);
		void* Symbolic, * Numeric;
		int result1 = umfpack_di_symbolic(3 * numT, 3 * numT, &col_ptr[0], &row_index[0], &entries[0], &Symbolic, NULL, NULL);
		int result2 = umfpack_di_numeric(&col_ptr[0], &row_index[0], &entries[0], Symbolic, &Numeric, NULL, NULL);


		std::vector<double> Qs(numT * 9);

		for (int c = 0; c < 3; c++) {
			double* b = &rhs_123[c][0];
			double* x = &Qs[c * 3 * numT];
			//int result3 = umfpack_di_solve(UMFPACK_A, &col_ptrQ[0], &row_indexQ[0], &entriesQ[0], x, b, NumericQ, NULL, NULL);
			int result3 = umfpack_di_solve(UMFPACK_A, &col_ptr[0], &row_index[0], &entries[0], x, b, Numeric, NULL, NULL);
		}

		for (UInt i1 = 0; i1 < numT; i1++) {
			double* r0 = &Qs[3 * i1];
			double* r1 = &Qs[3 * i1 + 3 * numT];
			double* r2 = &Qs[3 * i1 + 6 * numT];
			D[i1].setRowI(vec3d(r0), 0);
			D[i1].setRowI(vec3d(r1), 1);
			D[i1].setRowI(vec3d(r2), 2);
		}

		umfpack_di_free_numeric(&Numeric);
		umfpack_di_free_symbolic(&Symbolic);
	}

	/// <summary>
	/// Compute the final vertices of a new instance given preliminiary vertices where the resulting triangle mesh would not be consistent.
	/// 
	/// See equation (9) in [1].
	/// </summary>
	/// <param name="rots">
	///   The matrix "R_l[k] S_k Q_k" for every triangle k.
	///   R_l[k] is the joint rotation of the join associated with the triangle k.
	/// </param>
	/// <seealso cref="reconstructModel"/>
	/// <returns>
	///   Final vertices of the new instance, so that the triangle mesh given by the computed vertices
	///   and the ambient topology from the loaded SCAPE model is consistent.
	/// </returns>
	std::vector<vec3d> decodeDefGradMesh(std::vector<SMat3D>& rots) const {

		int numV = (int)default_mesh_vertices.size();

		// Build system matrix
#ifdef SAVE_MATRIX
		RandAccessCompRowMatrix<double> A(numV - 1, numV - 1);
#endif
		std::vector<double> rhs_x(numV);
		std::vector<double> rhs_y(numV);
		std::vector<double> rhs_z(numV);

		for (int i1 = 0; i1 < numT; i1++) {
			const int& v0 = default_mesh_triangles[i1].x;
			const int& v1 = default_mesh_triangles[i1].y;
			const int& v2 = default_mesh_triangles[i1].z;
			const SMat3D& R = rots[i1];
			vec3d vdash1 = default_mesh_vertices[v1] - default_mesh_vertices[v0];
			vec3d vdash2 = default_mesh_vertices[v2] - default_mesh_vertices[v0];
			vec3d Rv1 = R.vecTrans(vdash1);
			vec3d Rv2 = R.vecTrans(vdash2);


#ifdef SAVE_MATRIX
			if (v0 > 0) {
				A.add(v0 - 1, v0 - 1, 2);
				if (v1 > 0) A.add(v0 - 1, v1 - 1, -1);
				if (v2 > 0) A.add(v0 - 1, v2 - 1, -1);
			}
			if (v1 > 0) {
				A.add(v1 - 1, v1 - 1, 1);
				if (v0 > 0) A.add(v1 - 1, v0 - 1, -1);
			}
			if (v2 > 0) {
				A.add(v2 - 1, v2 - 1, 1);
				if (v0 > 0) A.add(v2 - 1, v0 - 1, -1);
			}
#endif

			rhs_x[v0] -= (Rv1.x + Rv2.x);
			rhs_y[v0] -= (Rv1.y + Rv2.y);
			rhs_z[v0] -= (Rv1.z + Rv2.z);
			rhs_x[v1] += Rv1.x;
			rhs_y[v1] += Rv1.y;
			rhs_z[v1] += Rv1.z;
			rhs_x[v2] += Rv2.x;
			rhs_y[v2] += Rv2.y;
			rhs_z[v2] += Rv2.z;
		}

#ifdef SAVE_MATRIX
		std::vector<double> entries;
		std::vector<int> row_index;
		std::vector<int> col_ptr;
		A.getMatrix(entries, row_index, col_ptr);
		void* Symbolic, * Numeric;
		int result1 = umfpack_di_symbolic(numV - 1, numV - 1, &col_ptr[0], &row_index[0], &entries[0], &Symbolic, NULL, NULL);
		int result2 = umfpack_di_numeric(&col_ptr[0], &row_index[0], &entries[0], Symbolic, &Numeric, NULL, NULL);
		umfpack_di_save_numeric(Numeric, "SCAPE_DGrad_numeric.bin");
		umfpack_di_save_symbolic(Symbolic, "SCAPE_DGrad_symbolic.bin");
		std::ofstream fout("matrixDGrad.bin", std::ios::binary);
		int numE = (int)entries.size();
		fout.write((char*)&numE, sizeof(int));
		fout.write((char*)&entries[0], sizeof(double) * numE);
		int numR = (int)row_index.size();
		fout.write((char*)&numR, sizeof(int));
		fout.write((char*)&row_index[0], sizeof(int) * numR);
		int numC = (int)col_ptr.size();
		fout.write((char*)&numC, sizeof(int));
		fout.write((char*)&col_ptr[0], sizeof(int) * numC);
		fout.close();
		umfpack_di_save_symbolic(Symbolic, "SCAPE_DGrad_symbolic.bin");
		umfpack_di_save_numeric(Numeric, "SCAPE_DGrad_numeric.bin");
#endif
		std::vector<double> res(3 * numV);


		//for (int i1 = 0; i1 < numV; i1++) std::cerr << rhs_x[i1] << std::endl;


		for (int c = 0; c < 3; c++) {
			double* b;
			if (c == 0) b = &rhs_x[1];
			else if (c == 1) b = &rhs_y[1];
			else b = &rhs_z[1];

			double* x = &res[c * numV + 1];
			int result3 = umfpack_di_solve(UMFPACK_A, &col_ptrDGrad[0], &row_indexDGrad[0], &entriesDGrad[0], x, b, NumericDGrad, NULL, NULL);
			//int result3 = umfpack_di_solve(UMFPACK_A, &col_ptr[0], &row_index[0], &entries[0], x, b, Numeric, NULL, NULL);
		}

#ifdef SAVE_MATRIX
		umfpack_di_free_numeric(&Numeric);
		umfpack_di_free_symbolic(&Symbolic);
#endif

		std::vector<vec3d> ret(numV);
#pragma omp parallel
		for (int i1 = 0; i1 < numV; i1++) {
			vec3d p(res[i1], res[numV + i1], res[2 * numV + i1]);
			ret[i1] = p;
		}
		return ret;
	}

	//private:

	void loadUMFPACKMatrices(const char* pathToMatrices) {

		int numE;
		int numR;
		int numC;
		//umfpack_di_load_numeric(&NumericQ, "SCAPE_Q_numeric.bin");
		//umfpack_di_load_symbolic(&SymbolicQ, "SCAPE_Q_symbolic.bin");
		//std::ifstream fin("matrixQ.bin", std::ios::binary);
		//fin.read((char*)&numE, sizeof(int));
		//entriesQ.resize(numE);
		//fin.read((char*)&entriesQ[0], sizeof(double)*numE);
		//fin.read((char*)&numR, sizeof(int));
		//row_indexQ.resize(numR);
		//fin.read((char*)&row_indexQ[0], sizeof(int)*numR);
		//fin.read((char*)&numC, sizeof(int));
		//col_ptrQ.resize(numC);
		//fin.read((char*)&col_ptrQ[0], sizeof(int)*numC);
		//fin.close();


		std::stringstream ss1;
		ss1 << pathToMatrices << "SCAPE_DGrad_numeric.bin";
		std::stringstream ss2;
		ss2 << pathToMatrices << "SCAPE_DGrad_symbolic.bin";
		std::stringstream ss3;
		ss3 << pathToMatrices << "matrixDGrad.bin";

		int res = umfpack_di_load_numeric(&NumericDGrad, const_cast<char*>(ss1.str().c_str()));
		if (res != 0) {
			std::cerr << "Error reading SCAPE_DGrad_numeric.bin" << std::endl;
			exit(1);
		}
		res = umfpack_di_load_symbolic(&SymbolicDGrad, const_cast<char*>(ss2.str().c_str()));
		if (res != 0) {
			std::cerr << "Error reading SCAPE_DGrad_symbolic.bin" << std::endl;
			exit(1);
		}
		std::ifstream fin2(ss3.str().c_str(), std::ios::binary);
		if (!fin2.good()) {
			std::cerr << "Error reading matrixDGrad.bin" << std::endl;
			exit(1);
		}
		fin2.read((char*)&numE, sizeof(int));
		entriesDGrad.resize(numE);
		fin2.read((char*)&entriesDGrad[0], sizeof(double) * numE);
		fin2.read((char*)&numR, sizeof(int));
		row_indexDGrad.resize(numR);
		fin2.read((char*)&row_indexDGrad[0], sizeof(int) * numR);
		fin2.read((char*)&numC, sizeof(int));
		col_ptrDGrad.resize(numC);
		fin2.read((char*)&col_ptrDGrad[0], sizeof(int) * numC);
		fin2.close();
	}

	void* SymbolicQ;
	void* NumericQ;
	std::vector<double> entriesQ;
	std::vector<int> row_indexQ;
	std::vector<int> col_ptrQ;

	void* SymbolicDGrad;
	void* NumericDGrad;
	std::vector<double> entriesDGrad;
	std::vector<int> row_indexDGrad;
	std::vector<int> col_ptrDGrad;


	//! Vertices of default mesh
	std::vector<vec3d> default_mesh_vertices;

	//! Triangles of default mesh
	std::vector<vec3i> default_mesh_triangles;

	//! The rigid parts of the model (stors for each part the triangles and the vertices that belong to the part)."
	Part rigidParts[numParts];

	//! The non-rigid pose deformation of a triangle depends on these two joints.
	std::vector< vec2i > triToJoints;

	//! The regressed linear functions for reconstucting the affine matrices Qi from the closest joint rotations.
	std::vector<perTriangleParams> perTriangleParameters;

	//! The PCA bases of the joint rotation vectors.
	std::vector< SMat3D > joint_PCA_bases;

	//! The joints of the skeleton.
	std::vector<Joint> joints;

	//! Specifies the part ID for each triangle
	std::vector<int> triIDtoPartID;

	//! The average deformation gradient for the pca basis.
	std::vector<double> avg_defgrad;

	//! Reconstruction the model from these SCAPE parameters will yield the average person in the average pose.
	std::vector<double> avg_SCAPE_coeffs;

	//! The pca basis for body shapes (Encoded as D matrices (deformation gradients)).
	PCA_Basis pcabasis;

	//! The pca basis for the joint angles.
	PCA_Basis pose_pcabasis;

	//! the average joint angles.
	std::vector<double> avg_jointangles;

	//! Lists for each triangle the 3 neighboring triangles.
	std::vector<vec3i> tri_neighbors;

	/// <summary>
	/// Triangles of of the topologically correct mesh. Use these when saving the mesh with <see cref="SimpleMesh"/>.
	/// This vector links three vertix indices together to one triangle in every element.
	/// </summary>
	/// 
	/// <remarks>
	/// All models (whether the default ['average'] one or any deformed derivative) consist of the same amount of vertices,
	/// which are on top of that also semantically identical.
	/// For example, if the vertex with index 1234 in the average model is the tip of the nose, then the vertex with index
	/// 1234 in any derived model will also reprsent the tip of the nose.
	/// This vector links three vertices together to one triangle in every element.
	/// In other words, every vec3i element contains three indices i1, i2, i3, so that the vertices with indices i1, i2, i3
	/// form a triangle.
	/// </remarks>
	/// 
	/// <see cref="SimpleMesh"/>
	std::vector<vec3i> topologically_correct_mesh_triangles;

};
