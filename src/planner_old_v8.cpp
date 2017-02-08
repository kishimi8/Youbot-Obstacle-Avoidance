
#include "ros/ros.h"
#include "nav_msgs/Odometry.h"
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <math.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Segment_2<K> Segment_2;
typedef CGAL::Line_2<K> Line_2;
typedef CGAL::Vector_2<K> Vector_2;


double crpdt(Point_2 a, Point_2 b);
double absolute(Point_2 a);
double dst_pt_lineseg(Segment_2 l, Point_2 pt, Point_2 &pt_cl);
double rep_potential_field(double dist_rep);
double att_potential_field(double dist_att);

class Listener
{
public:
 double x;
 double y;
 double thx; double thy; double thz; double thw;
 double x_fl; double y_fl; 
 double x_bl; double y_bl; 
 double x_fr; double y_fr; 
 double x_br; double y_br; 
Listener()

	{x=0; y=0; thx=0; thy=0; thz = 0; thw =0; 
	x_fr = x + 0.29; // Front right Corner
y_fr = y - 0.19;

x_br = x - 0.29; // Back right Corner
y_br = y - 0.19;

x_fl = x + 0.29; // Front left Corner
y_fl = y + 0.19;

x_bl = x - 0.29; // Back left Corner
y_bl = y + 0.19; }
void chatterCallback(const nav_msgs::Odometry::ConstPtr& msg)
{
//ROS_INFO("Seq: [%d]", msg->header.seq);
/*  ROS_INFO("Position-> x: [%f], y: [%f], z: [%f]", msg->pose.pose.position.x,msg->pose.pose.position.y, msg->pose.pose.position.z);
  ROS_INFO("Orientation-> x: [%f], y: [%f], z: [%f], w: [%f]", msg->pose.pose.orientation.x, msg->pose.pose.orientation.y, msg->pose.pose.orientation.z, msg->pose.pose.orientation.w);
  ROS_INFO("Vel-> Linear: [%f], Angular: [%f]", msg->twist.twist.linear.x,msg->twist.twist.angular.z);
*/
x = msg->pose.pose.position.x;
y = msg->pose.pose.position.x;
thx = msg->pose.pose.orientation.x;
thy = msg->pose.pose.orientation.y;
thz = msg->pose.pose.orientation.z;
thw = msg->pose.pose.orientation.w;
// Control points for the robot
// We have chosen the corner points.. Remember the dimensions of youbot are:
// Total Length: 580 mm, Total Width: 380 mm

x_fr = x + 0.29; // Front right Corner
y_fr = y - 0.19;

x_br = x - 0.29; // Back right Corner
y_br = y - 0.19;

x_fl = x + 0.29; // Front left Corner
y_fl = y + 0.19;

x_bl = x - 0.29; // Back left Corner
y_bl = y + 0.19;
}

};


// Projection on a line and a line segment

double dst_pt_lineseg (Segment_2 s, Point_2 pt, Point_2 &pt_cl)

{
double d1, d2, d, distance;
Point_2 proj_pt;
Line_2 l = s.supporting_line();
proj_pt = l.projection(pt); 


if (s.has_on(proj_pt))
{
//std::cout<<proj_pt.x()<<"\t"<<proj_pt.y()<<std::endl;
pt_cl = proj_pt;
distance = (double) sqrt(squared_distance (proj_pt, pt));
}
else
{
d1 = squared_distance(proj_pt,s.vertex(0));
d2 = squared_distance(proj_pt,s.vertex(1));
d = std::min(d1,d2);

	if (d == d1){
//		std::cout<<s.vertex(0)<<std::endl;
		distance = (double) sqrt(squared_distance(s.vertex(0),pt));
		pt_cl = s.vertex(0);
//		std::cout<<distance<<std::endl;
		    }
	else
		{
	//	std::cout<<s.vertex(1)<<std::endl;
		pt_cl = s.vertex(1);
		distance = (double) sqrt(squared_distance(s.vertex(1),pt));
		}
	
}
//std::cout<<distance<<std::endl;
return distance;
}

double dst_pt_pgn(Polygon_2 pgn, Point_2 pt, Point_2 &pt_cl)
{
double distance;
double distance1, distance2, distance3, distance4;
double min1, min2;
Point_2 p1, p2, p3, p4;

Segment_2 s1, s2, s3, s4;
s1 = pgn.edge(0);
s2 = pgn.edge(1);s3 = pgn.edge(2);s4 = pgn.edge(3);
distance1 = dst_pt_lineseg(s1,pt,p1);
distance2 = dst_pt_lineseg(s2,pt,p2);
distance3 = dst_pt_lineseg(s3,pt,p3);
distance4 = dst_pt_lineseg(s4,pt,p4);
min1 = std::min(distance1, distance2);
min2 = std::min(distance3, distance4);
distance = std::min(min1,min2);
if (distance == min1)
{
	if (min1 == distance1)
		pt_cl = p1;
	else
		pt_cl = p2;
}
else
{
	if (min2 == distance3)
		pt_cl = p3;
	else
		pt_cl = p4;
}

//std::cout<<distance1<<"\t"<<distance2<<"\t"<<distance3<<"\t"<<distance4<<std::endl;

return distance;
}


double dst_from_obs(Point_2 p_ctrl, Polygon_2 pgn_1, Polygon_2 pgn_2, Polygon_2 pgn_3, Polygon_2 pgn_4, Point_2 &pt_cl)
{
// Gives distance to the nearest obstacle as well as the position of the nearest point on the obstacle
double distance;
double d1, d2, d3, d4;
double min1, min2;
Point_2 p1, p2, p3, p4;

//Point_2 fl(p.x_fl, p.y_fl);
//Point_2 fr(p.x_fr, p.y_fr);
//Point_2 bl(p.x_bl, p.y_bl);
//Point_2 br(p.x_br, p.y_br);

d1 = dst_pt_pgn(pgn_1, p_ctrl, p1);
d2 = dst_pt_pgn(pgn_2, p_ctrl, p2);
d3 = dst_pt_pgn(pgn_3, p_ctrl, p3);
d4 = dst_pt_pgn(pgn_4, p_ctrl, p4);

min1 = std::min(d1,d2);
min2 = std::min(d3,d4);

// Getting Minimum distance of all the obstacles
distance = std::min(min1,min2);

// Getting Closest Point of all the obstacles
if (distance == min1)
{
	if (min1 == d1)
		pt_cl = p1;
	else
		pt_cl = p2;
}
else
{
	if (min2 == d3)
		pt_cl = p3;
	else
		pt_cl = p4;
}


std::cout<<distance<<std::endl;
std::cout<<pt_cl<<std::endl;
return distance;
}

void planner(Listener &p, Polygon_2 pgn_1, Polygon_2 pgn_2, Polygon_2 pgn_3, Polygon_2 pgn_4, Point_2 goal)
{
  double dist_rep;
  double att_gain = 2;
  double att_thresh = 0.5;
  double dist_att_fl,dist_att_fr,dist_att_bl,dist_att_br;
  double dist_rep_fl,dist_rep_fr,dist_rep_bl,dist_rep_br;
  Point_2 fl(p.x_fl, p.y_fl);
  Point_2 fr(p.x_fr, p.y_fr);
  Point_2 bl(p.x_bl, p.y_bl);
  Point_2 br(p.x_br, p.y_br);
  Point_2 nearest; 

  dist_rep_fl = dst_from_obs(fl, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
  dist_rep_fr = dst_from_obs(fr, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
  dist_rep_bl = dst_from_obs(bl, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
  dist_rep_br = dst_from_obs(br, pgn_1,pgn_2,pgn_3,pgn_4,nearest);


  dist_att_fl = (double) sqrt(squared_distance(fl, goal));
  dist_att_fr = (double) sqrt(squared_distance(fr, goal));
  dist_att_bl = (double) sqrt(squared_distance(bl, goal));
  dist_att_br = (double) sqrt(squared_distance(br, goal));

 // Attractive Potential Field  

//  dU_att_fl = 
//  dU_att_fr = 
//  dU_att_bl = 
//  dU_att_br = 

 // Repulsive Potential Field

 // dU_rep_fl = 
 // dU_rep_fr = 
 // dU_rep_bl = 
 // dU_rep_br = 

 // Forces in Workspace
 // dU_att = sum (dU_att_fl)+sum (dU_att_fr)+sum (dU_att_bl)+sum (dU_att_br);
 // dU_rep = 
 // Forces in ConfigSpace

//    u = config_force(dU_att_fl, p.jacobian_fl) + config_force(dU_att_fr, p.jacobian_fr) + config_force(dU_att_bl, p.jacobian_bl) + config_force(dU_att_br, p.jacobian_br);
//   u has three components here
}



double att_potential_field(double dist_att)
{
  double dist_rep;
  double att_gain = 2;
  double att_thresh = 0.5;
  double dU_att, U_att;
  if (dist_att < att_thresh) 
		{
   U_att =  0.5 * att_gain * pow(dist_att,2);
   dU_att = att_gain * dist_att;
 		}
  else
		{
   U_att = att_thresh * att_gain * dist_att - 0.5 * att_gain * pow(att_thresh,2);
   dU_att = (double) (att_thresh * att_gain * dist_att)/dist_att;
		}

return dU_att;

}

double rep_potential_field(double dist_rep)
{
  double rep_gain = 2;
  double rep_thresh = 0.5;
  double dU_rep;

  if (dist_rep < rep_thresh) 
		{

 		}
  else
		{

		}
return dU_rep;
}

int main(int argc, char** argv)
{
  double dist_rep, dist_att;
  Point_2 goal(5,5);
  geometry_msgs::Twist vel;
  Listener pose;
  ros::init(argc, argv, "subscriber_node"); 
  ros::NodeHandle n1; 
  ros::init(argc, argv, "publisher_node");  
  ros::NodeHandle n2; 
  ros::Publisher pub = n1.advertise<geometry_msgs::Twist>("cmd_vel",1000);
  ros::Subscriber sub = n1.subscribe("odom",1000,&Listener::chatterCallback,&pose); 
  ros::Rate loop_rate(5);


/* 
Need to save the obstacles in the form of polygons here
*/
   Point_2 o11(2,2), o12(3,2), o13(3,3), o14(2,3); // obs_1 red in color
   Point_2 obs_1[] = {o11,o12,o13,o14};
   Polygon_2 pgn_1(obs_1, obs_1+4);

   Point_2 o21(4.625,0.625), o22(5.375,0.625), o23(5.375,1.375), o24(4.625,1.375); // obs_2 green in color
   Point_2 obs_2[] = {o21,o22,o23,o24};
   Polygon_2 pgn_2(obs_2, obs_2+4);

   Point_2 o31(2.75,4.75), o32(3.25,4.75), o33(3.25,5.25), o34(2.75,5.25); // obs_3 blue in color
   Point_2 obs_3[] = {o31,o32,o33,o34};
   Polygon_2 pgn_3(obs_3, obs_3+4);

   Point_2 o41(-0.125,2.875), o42(0.125,2.875), o43(0.125,3.125), o44(-0.125,3.125); // obs_4 purple in color
   Point_2 obs_4[] = {o41,o42,o43,o44};
   Polygon_2 pgn_4(obs_4, obs_4+4);



   Point_2 nearest;
while(ros::ok())
 { 
  ros::spinOnce();

  pub.publish(vel);

//Point_2 m(pose.x_fl,pose.y_fl);
//Point_2 m(5,5);
//Point_2 res;
//Segment_2 s(o11,o12);
//distance = dst_pt_lineseg(s, m);
//distance = dst_pt_pgn(pgn_1, m, res);

planner(pose,pgn_1,pgn_2,pgn_3,pgn_4, goal);
//std::cout<<distance<<"\t"<<res<<std::endl;

loop_rate.sleep();
}
  
  vel.linear.x = 0.0;
  vel.linear.y = 0.0;
  vel.linear.z = 0;

 

  return 0;
}



