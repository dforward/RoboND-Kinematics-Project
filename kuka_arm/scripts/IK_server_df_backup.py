#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:

        ### Your FK code here

	# Create symbols
   	# Twist Angles
    	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
    	# Link Lengths
    	a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    	# Link Offsets
    	d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')    				
    	# Joint Angles
    	q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

	
	# Create Modified DH parameters
	# DH Parameter Table
    	dh_param_table = {alpha0: 0,     a0: 0,      d1: 0.75,  q1: q1,
	              alpha1: -pi/2, a1: 0.35,   d2: 0,     q2: -pi/2 + q2,
		      alpha2: 0,     a2: 1.25,   d3: 0,     q3: q3,
	              alpha3: -pi/2, a3: -0.054, d4: 1.5,   q4: q4,
		      alpha4: pi/2,  a4: 0,      d5: 0,     q5: q5,
		      alpha5: -pi/2, a5: 0,      d6: 0,     q6: q6,
		      alpha6: 0,     a6: 0,      d7: 0.303, q7: 0}	


	# Define Modified DH Transformation matrix
	def dh_trans_matrix(alpha, a, d, q):
		dh_tm = Matrix([[cos(q),            -sin(q),           0,           a            ],
     			[sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
			[sin(q)*sin(alpha), cos(q)*sin(alpha), cos(alpha),  cos(alpha)*d ],
			[0,                 0,                 0,           1            ]])
		return dh_tm
	
	# Create individual transformation matrices
	t01 = dh_trans_matrix(alpha0, a0, d1, q1).subs(dh_param_table)
    	t12 = dh_trans_matrix(alpha1, a1, d2, q2).subs(dh_param_table)
    	t23 = dh_trans_matrix(alpha2, a2, d3, q3).subs(dh_param_table)
    	t34 = dh_trans_matrix(alpha3, a3, d4, q4).subs(dh_param_table)
    	t45 = dh_trans_matrix(alpha4, a4, d5, q5).subs(dh_param_table)
    	t56 = dh_trans_matrix(alpha5, a5, d6, q6).subs(dh_param_table)
    	t6ee = dh_trans_matrix(alpha6, a6, d7, q7).subs(dh_param_table)
	
	# Homogeneous Transform Matrix from base link to End Effector
	t0ee = t01 * t12 * t23 * t34 * t45 * t56 * t6ee

	# Corrections to account for difference between gripper in URDF and DH Convention	

    	R_y = Matrix([[ cos(-pi/2),        0,  sin(-pi/2)],
                      [       0,        1,        0],
                      [-sin(-pi/2),        0,  cos(-pi/2)]])

    	R_z = Matrix([[ cos(pi), -sin(pi),        0],
                      [ sin(pi),  cos(pi),        0],
                      [ 0,              0,      1]])

    	R_ee = R_z * R_y

    	# Correct for End Effector rotation error
    	t0ee_corrected = t0ee * R_ee
    
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here

	    # Get Gripper/End Effector Rotation Matrix
    	    r, p, y = symbols('r, p, y')

    	    roll_rot = Matrix([[ 1,              0,        0],
    	                  [ 0,        cos(r), -sin(r)],
              	          [ 0,        sin(r),  cos(r)]])

    	    pitch_rot = Matrix([[ cos(p),        0,  sin(p)],
                          [       0,        1,        0],
                          [-sin(p),        0,  cos(p)]])

    	    yaw_rot = Matrix([[ cos(y), -sin(y),        0],
                          [ sin(y),  cos(y),        0],
                          [ 0,              0,      1]])

    	    R_ee = yaw_rot * pitch_rot * roll_rot

    	    # Correct for End Effector rotation error
    	    R_error = yaw_rot.subs(y, radians(180)) * pitch_rot.subs(p, radians(-90))
    	    R_ee = R_ee * R_error	

	    R_ee = R_ee.subs({'r': roll, 'p': pitch, 'y': yaw})

	    EE_pos = Matrix([[px],
     	    [py],
	    [pz]])

	    wrist_center = EE_pos - (0.303) * R_ee[:,2]

	    # Calculate joint angles using Geometric IK method
	   
	    # Calculate Theta1 using the x and y positions of the wrist center	
    	    theta1 = atan2(wrist_center[1], wrist_center[0])

	    # Use Law of Cosines to solve triangle and get theta2 and theta3
    	    A = 1.5
    	    B = sqrt(pow((sqrt(wrist_center[0] * wrist_center[0] + wrist_center[1] * wrist_center[1]) - 0.35), 2) + pow((wrist_center[2] - 0.75), 2))
            C = 1.25    
    	    a_angle = acos((B*B + C*C - A*A) / (2*B*C))
    	    b_angle = acos((A*A + C*C - B*B) / (2*A*C))
    	    c_angle = acos((A*A + B*B - C*C) / (2*A*B))

	    # Calculate theta2 and theta3
    	    theta2 = pi/2 - a_angle - atan2(wrist_center[2] - 0.75, sqrt(wrist_center[0]*wrist_center[0] + wrist_center[1]*wrist_center[1]) - 0.35)
     	    theta3 = pi/2 - (b_angle + 0.036)

	    # Rotation Matrix from base link to third link
    	    R03 = t01[0:3,0:3] * t12[0:3,0:3] * t23[0:3,0:3]
    	    R03 = R03.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
    	    # Get rotation matrix from 3 to 6
    	    R36 = R03.inv("LU") * R_ee

    	    # Get the Eurler angles from the rotation Matrix
    	    theta4 = atan2(R36[2,2], -R36[0,2])
    	    theta5 = atan2(sqrt(R36[0,2] * R36[0,2] + R36[2,2] * R36[2,2]), R36[1,2])
    	    theta6 = atan2(-R36[1,1], R36[1,0])

            ###

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
