using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

/// <summary>
/// Basic constraint which needs to be dropped on a sphere.
/// We only implement spherical joints. Other constraints would require a generalization of the constraint interface.
/// </summary>
public class Constraint : MonoBehaviour
{
	/// <summary>
	/// Default constructor. All zero. 
	/// </summary>
	public Constraint ()
	{
		this.m_manager = null;
	}

	#region EditorVariables

	public RigidBody BodyA;
	public RigidBody BodyB;

	#endregion

	#region OtherVariables

	private PhysicsManager m_manager;
	private int m_index;
	//Index of the constraint into the global vector of constraints
	private Vector3 m_pA;
	//Constraint point in the local reference frame of bodyA
	private Vector3 m_pB;
	//Constraint point in the local reference frame of bodyB
	// If BodyA or BodyB is not defined, m_pA or m_pB stores the global coordinates of the constraint point

	#endregion

	#region MonoBehaviour

	// Nothing to do here

	#endregion

	#region OtherMethods

	public void initialize (int index, PhysicsManager manager)
	{
		m_manager = manager;

		m_index = index;

		// Get the center of the sphere as the constraint point and transform it to the local frames of the bodies
		Transform xform = this.GetComponent<Transform> ();

		if (xform != null) {
			Vector3 p = xform.localPosition;

			if (BodyA != null) {
				m_pA = BodyA.PointGlobalToLocal (p);
			} else {
				m_pA = p;
			}

			if (BodyB != null) {
				m_pB = BodyB.PointGlobalToLocal (p);
			} else {
				m_pB = p;
			}
		}
	}

	public int getSize ()
	{
		return 3;
	}

	public void updateScene ()
	{
		// Apply the average position to the mesh
		this.GetComponent<Transform> ().position =
            0.5f * ((BodyA ? BodyA.PointLocalToGlobal (m_pA) : m_pA) + (BodyB ? BodyB.PointLocalToGlobal (m_pB) : m_pB));
	}

	public void addForces ()
	{
		if (BodyA == null && BodyB == null)
			return;
			
		// Compute constraint points global coordinates
		Vector3 globalMPA = getGlobalPosition (BodyA, m_pA);
		Vector3 globalMPB = getGlobalPosition (BodyB, m_pB);

		// Constraint
		float C = (globalMPA - globalMPB).magnitude;

		// Add forces to each rigid body, if present
		if (BodyA != null)
			BodyA.addSoftConstraintForces (-m_manager.K * C * (globalMPA - globalMPB).normalized, m_pA);

		if (BodyB != null)
			BodyB.addSoftConstraintForces (-m_manager.K * C * (globalMPB - globalMPA).normalized, m_pB);

	}

	// Computes the constraint contact position in global coordinates
	private Vector3 getGlobalPosition (RigidBody body, Vector3 constraintPoint)
	{
		if (body == null)
			return transform.position;
		else
			return body.PointLocalToGlobal (constraintPoint);
	}

	public void getConstraintVector (VectorXD C0)
	{
		if (BodyA == null && BodyB == null)
			return;
		
		Vector3 globalA = getGlobalPosition (BodyA, m_pA);
		Vector3 globalB = getGlobalPosition (BodyB, m_pB);

		C0.SetSubVector (m_index, 3, Utils.ToVectorXD (globalA - globalB));
	}

	public void getConstraintJacobian (MatrixXD J)
	{
		getConstraintJacobianForBody (BodyA, m_pA, J, 1.0f);
		getConstraintJacobianForBody (BodyB, m_pB, J, -1.0f);
	}

	// Compute the velocity and angular velocity derivatives for a body, if exists
	private void getConstraintJacobianForBody (RigidBody body, Vector3 localPoint, MatrixXD J, float sign)
	{
		if (body) {
			int vi = body.getSimIndex ();
			int wi = vi + 3;

			MatrixXD Jv = DenseMatrixXD.CreateIdentity (3);
			MatrixXD Jw = Utils.Skew (-body.VecLocalToGlobal (localPoint));

			J.SetSubMatrix (m_index, 3, vi, 3, (Jv * sign));
			J.SetSubMatrix (m_index, 3, wi, 3, (Jw * sign));
		} 
	}

	#endregion
}
