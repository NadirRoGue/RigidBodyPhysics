using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;


/// <summary>
/// Basic physics manager capable of simulating a given ISimulable
/// implementation using diverse integration methods: explicit,
/// implicit, Verlet and semi-implicit.
/// </summary>
public class PhysicsManager : MonoBehaviour
{
	/// <summary>
	/// Default constructor. Zero all. 
	/// </summary>
	public PhysicsManager ()
	{
		this.Paused = true;
		this.OneStep = false;
		this.TimeStep = 0.01f;
		this.Gravity = new Vector3 (0.0f, -9.81f, 0.0f);
		this.MethodConstraints = ConstraintMethod.Hard;
	}

	/// <summary>
	/// Constraint method.
	/// </summary>
	public enum ConstraintMethod
	{
		Soft = 0,
		Hard = 1}

	;

	#region InEditorVariables

	public bool Paused;
	public bool OneStep;
	public float TimeStep;
	public Vector3 Gravity;
	public float K;
	public float Damping;
	public List<GameObject> SimulableObjects;
	public List<GameObject> Constraints;
	public ConstraintMethod MethodConstraints;

	#endregion

	#region OtherVariables

	private List<ISimulable> m_objs;
	private List<Constraint> m_constraints;
	private int m_numdofs;
	private int m_numcs;

	#endregion

	#region MonoBehaviour

	public void Start ()
	{
		int index = 0;

		m_objs = new List<ISimulable> (SimulableObjects.Capacity);

		foreach (GameObject obj in SimulableObjects) {
			ISimulable simobj = obj.GetComponent<ISimulable> ();

			if (simobj != null) {
				m_objs.Add (simobj);

				// Initialize simulable model
				simobj.initialize (index, this);

				// Retrieve pos and vel size
				index += simobj.getNumDof ();
			}
		}

		m_numdofs = index;

		index = 0;

		m_constraints = new List<Constraint> (Constraints.Capacity);

		foreach (GameObject obj in Constraints) {
			Constraint constraint = obj.GetComponent<Constraint> ();

			if (constraint != null) {
				m_constraints.Add (constraint);

				// Initialize constraint
				constraint.initialize (index, this);

				// Retrieve constraint size
				index += constraint.getSize ();
			}
		}

		m_numcs = index;
	}

	public void Update ()
	{
		if (Input.GetKeyUp (KeyCode.O))
			this.OneStep = !this.OneStep;

		if (Input.GetKeyUp (KeyCode.P))
			this.Paused = !this.Paused;
	}

	public void FixedUpdate ()
	{
		if (this.Paused && !this.OneStep)
			return; // Not simulating

		if (this.OneStep) // One!
			this.OneStep = false;

		// Integration method
		switch (MethodConstraints) {
		case ConstraintMethod.Hard:
			HardConstraintsStep ();
			break;
		case ConstraintMethod.Soft:
		default:
			SoftConstraintsStep ();
			break;
		}

		// Update visual elements
		foreach (ISimulable obj in m_objs) {
			obj.updateScene ();
		}

		foreach (Constraint constraint in m_constraints) {
			constraint.updateScene ();
		}
	}

	#endregion

	#region OtherMethods

	private void HardConstraintsStep ()
	{
		// Vectors and matrices to store data
		VectorXD C0 = new DenseVectorXD (m_numcs);
		VectorXD f = new DenseVectorXD (m_numdofs);
		VectorXD v = new DenseVectorXD (m_numdofs);
		MatrixXD M = new DenseMatrixXD (m_numdofs, m_numdofs);
		MatrixXD J = new DenseMatrixXD (m_numcs, m_numdofs);

		MatrixXD System = new DenseMatrixXD (m_numcs + m_numdofs, m_numcs + m_numdofs);
		VectorXD Independent = new DenseVectorXD (m_numdofs + m_numcs);

		// Compute forces and retrieve the rigid bodies data
		foreach (ISimulable i in m_objs) {
			
			i.clearForcesAndMatrices ();
			i.addForcesAndMatrices ();

			i.getForceVector (f);
			i.getVelocityVector (v);
			i.getMassMatrix (M);
		}

		// Compute and retrieve restrictions and Jacobians
		foreach (Constraint c in m_constraints) {
			c.getConstraintVector (C0);
			c.getConstraintJacobian (J);
		}

		// Set up left-side system
		System.SetSubMatrix (0, m_numdofs, 0, m_numdofs, M);
		System.SetSubMatrix (m_numdofs, m_numcs, 0, m_numdofs, J);
		System.SetSubMatrix (0, m_numdofs, m_numdofs, m_numcs, J.Transpose ());
		System.SetSubMatrix (m_numdofs, m_numcs, m_numdofs, m_numcs, DenseMatrixXD.Create (m_numcs, m_numcs, 0.0));

		// Set up right-side system
		VectorXD b = M * v + TimeStep * f;
		VectorXD AtC0 = (-1.0f / TimeStep) * C0;
		Independent.SetSubVector (0, m_numdofs, b);
		Independent.SetSubVector (m_numdofs, m_numcs, AtC0);

		// Solve system
		VectorXD newVelocities = System.Solve (Independent);

		// Update bodies
		foreach (ISimulable i in m_objs) {
			i.setVelocityVector (newVelocities);
			i.advancePosition ();
		}
	}

	private void SoftConstraintsStep ()
	{
		// Clear previous data
		foreach (ISimulable i in m_objs) {
			i.clearForcesAndMatrices ();
		}

		// Add forces caused by constraints
		foreach (Constraint c in m_constraints) {
			// Compute forces and torque caused by soft constraints
			c.addForces ();
		}
			
		// Integrate rigid body
		foreach (ISimulable i in m_objs) {
			RigidBody rb = i as RigidBody;

			// Compute forces and torque caused by gravity and damping
			i.addForcesAndMatrices ();

			// Needed data
			Vector3 F = rb.getForce ();
			Vector3 T = rb.getTorque ();
			Vector3 v = rb.getVelocity ();
			Vector3 w = rb.getAngularVelocity ();
			MatrixXD M = rb.getInertiaMatrix ();

			// Integrate velocity
			Vector3 newV = v + TimeStep * F / rb.mass;
			rb.setVelocity (newV);
				
			// Integrate angular velocity
			VectorXD L = M * Utils.ToVectorXD (w);
			L += Utils.ToVectorXD (TimeStep * T);
			VectorXD newW = M.Inverse () * L;
			rb.setAngularVelocity (Utils.ToVector3 (newW));

			// Update position and rotation
			i.advancePosition ();
		}
	}

	#endregion

}
