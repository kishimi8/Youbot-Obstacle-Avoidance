#include <CGAL/Random.h>
#include "ros/ros.h"
#include "nav_msgs/Odometry.h"
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <math.h>
#include <CGAL/Cartesian_matrix.h>
#include <iterator>
#include "cmatrix"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Segment_2<K> Segment_2;
typedef CGAL::Line_2<K> Line_2;
typedef CGAL::Vector_2<K> Vector_2;
typedef CGAL::Vector_3<K> Vector_3;
typedef techsoft::matrix<double> Matrix;

double crpdt(Point_2 a, Point_2 b);
double absolute(Point_2 a);
double dst_pt_lineseg(Segment_2 l, Point_2 pt, Point_2 &pt_cl);
double rep_potential_field(double dist_rep);
//Matrix att_potential_field(double dist_att, Point_2 goal, Point_2 ctrl_pt);
Point_2 att_potential_field(double dist_att, Point_2 goal, Point_2 ctrl_pt);
//Matrix config_force(Matrix wksp_force, Matrix jac);
void config_force(Matrix wksp_force, Matrix jac, Matrix &u);

class Listener
{
public:
 double x;
 double y;
 double thx; double thy; double thz; double thw; double th;
 double x_fl; double y_fl; 
 double x_bl; double y_bl; 
 double x_fr; double y_fr; 
 double x_br; double y_br; 
 ~Listener();
Listener()

	{x=0; y=0; thx=0; thy=0; thz = 0; thw =0;  th=0;
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


x = msg->pose.pose.position.x;
y = msg->pose.pose.position.x;
th = msg->twist.twist.angular.z;
// Control points for the robot
// We have chosen the corner points.. Remember the dimensions of youbot are:
// Total Length: 580 mm, Total Width: 380 mm

x_fr = x -0.29*sin(th)+0.19*cos(th); // Front right Corner
y_fr = y +0.29*cos(th)+0.19*sin(th);;

x_br = x +0.29*sin(th)+0.19*cos(th); // Back right Corner
y_br = y -0.29*cos(th)+0.19*sin(th);

x_fl = x -0.29*sin(th)-0.19*cos(th); // Front left Corner
y_fl = y +0.29*cos(th)-0.19*sin(th);

x_bl = x +0.29*sin(th)-0.19*cos(th); // Back left Corner
y_bl = y -0.29*cos(th)-0.19*sin(th);
}

Point_2 jacobian_fl()
{

std::cout<<"came to jac"<<std::endl;
double r11, r21;
double j1, j2;
r11 = 0.29;
r21 = 0.19;

j1 = -r11*sin(th)-r21*cos(th);
j2 = r11*cos(th)-r21*sin(th);

Point_2 result(j1,j2);

return result;
}

Point_2 jacobian_fr()
{

std::cout<<"came to jac"<<std::endl;
double r11, r21;
double j1, j2;
r11 = 0.29;
r21 = -0.19;

j1 = -r11*sin(th)-r21*cos(th);
j2 = r11*cos(th)-r21*sin(th);

Point_2 result(j1,j2);

return result;


}

Point_2 jacobian_bl()
{
std::cout<<"came to jac"<<std::endl;
double r11, r21;
double j1, j2;
r11 = -0.29;
r21 = 0.19;


j1 = -r11*sin(th)-r21*cos(th);
j2 = r11*cos(th)-r21*sin(th);

Point_2 result(j1,j2);

return result;
}

Point_2 jacobian_br()
{
std::cout<<"came to jac"<<std::endl;
double r11, r21;
double j1, j2;
r11= -0.29;
r21 = -0.19;


j1 = -r11*sin(th)-r21*cos(th);
j2 = r11*cos(th)-r21*sin(th);

Point_2 result(j1,j2);

return result;

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

pt_cl = proj_pt;
distance = (double) sqrt(squared_distance (proj_pt, pt));
}
else
{
d1 = squared_distance(proj_pt,s.vertex(0));
d2 = squared_distance(proj_pt,s.vertex(1));
d = std::min(d1,d2);

	if (d == d1){

		distance = (double) sqrt(squared_distance(s.vertex(0),pt));
		pt_cl = s.vertex(0);

		    }
	else
		{

		pt_cl = s.vertex(1);
		distance = (double) sqrt(squared_distance(s.vertex(1),pt));
		}
	
}

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


return distance;
}


double dst_from_obs(Point_2 p_ctrl, Polygon_2 pgn_1, Polygon_2 pgn_2, Polygon_2 pgn_3, Polygon_2 pgn_4, Point_2 &pt_cl)
{
// Gives distance to the nearest obstacle as well as the position of the nearest point on the obstacle
double distance;
double d1, d2, d3, d4;
double min1, min2;
Point_2 p1, p2, p3, p4;

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


//std::cout<<distance<<std::endl;
//std::cout<<pt_cl<<std::endl;
return distance;
}

Vector_3 planner(Listener &p, Polygon_2 pgn_1, Polygon_2 pgn_2, Polygon_2 pgn_3, Polygon_2 pgn_4, Point_2 goal)
{

  double dist_rep;
  double att_gain = 2;
  double att_thresh = 0.5;
  double dist_att_fl,dist_att_fr,dist_att_bl,dist_att_br;
  double dist_rep_fl,dist_rep_fr,dist_rep_bl,dist_rep_br;
  
  //Matrix j1(2,3,1.0); Matrix j2(2,3,1.0); Matrix j3(2,3,1.0); Matrix j4(2,3,1.0);
/*  Matrix u1(3,1,1.0); Matrix u2(3,1,1.0); Matrix u3(3,1,1.0); Matrix u4(3,1,1.0);
  Matrix dU_att_fl(2,1,1.0); Matrix dU_att_fr(2,1,1.0); Matrix dU_att_bl(2,1,1.0); Matrix dU_att_br(2,1,1.0);
  Matrix dU_att(2,1,1.0);
  Matrix u_att(3,1,1.0), u_rep(3,1,1.0);
  */
/*    Matrix j_fl(2,3,1.0); Matrix j_fr(2,3,1.0); Matrix j_bl(2,3,1.0); Matrix j_br(2,3,1.0);
   j_fl[1][2] = 0; j_fl[2][1] = 0;
   j_fr[1][2] = 0; j_fr[2][1] = 0;
   j_bl[1][2] = 0; j_bl[2][1] = 0;
   j_br[1][2] = 0; j_br[2][1] = 0;
*/
  Point_2 fl(p.x_fl, p.y_fl);
  Point_2 fr(p.x_fr, p.y_fr);
  Point_2 bl(p.x_bl, p.y_bl);
  Point_2 br(p.x_br, p.y_br);
  Point_2 nearest; 

///  dist_rep_fl = dst_from_obs(fl, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
//  dist_rep_fr = dst_from_obs(fr, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
//  dist_rep_bl = dst_from_obs(bl, pgn_1,pgn_2,pgn_3,pgn_4,nearest);
//  dist_rep_br = dst_from_obs(br, pgn_1,pgn_2,pgn_3,pgn_4,nearest);


  dist_att_fl = (double) sqrt(squared_distance(fl, goal));
  dist_att_fr = (double) sqrt(squared_distance(fr, goal));
  dist_att_bl = (double) sqrt(squared_distance(bl, goal));
  dist_att_br = (double) sqrt(squared_distance(br, goal));
 
 // Attractive Potential Field  
//  OR
 // Forces in Workspace


    Point_2 dU_att_fl, dU_att_fr, dU_att_bl, dU_att_br;
    dU_att_fl = att_potential_field(dist_att_fl,goal, fl);
    dU_att_fr = att_potential_field(dist_att_fr,goal, fr);
    dU_att_bl = att_potential_field(dist_att_bl,goal, bl);
    dU_att_br = att_potential_field(dist_att_br,goal, br);
    /*Matrix w_fl(2,1,1.0);Matrix w_fr(2,1,1.0);   Matrix w_bl(2,1,1.0);    Matrix w_br(2,1,1.0);

    w_fl[1][1] = dU_att_fl.x();    w_fl[2][1] = dU_att_fl.y();
    w_fl[1][1] = dU_att_fr.x();    w_fl[2][1] = dU_att_fr.y();
    w_fl[1][1] = dU_att_bl.x();    w_fl[2][1] = dU_att_bl.y();
    w_fl[1][1] = dU_att_br.x();    w_fl[2][1] = dU_att_br.y();

*/
 // Repulsive Potential Field

 // dU_rep_fl = 
 // dU_rep_fr = 
 // dU_rep_bl = 
 // dU_rep_br = 

 // Forces in ConfigSpace

    Point_2 j1, j2, j3, j4;     
    j1 = p.jacobian_fl();
    j2 = p.jacobian_fr();
    j3 = p.jacobian_bl();
    j4 = p.jacobian_br();
    double u1, u2, u3;
    u1 = dU_att_fl.x() + dU_att_fr.x() + dU_att_bl.x() + dU_att_br.x();
    u2 = dU_att_fl.y() + dU_att_fr.y() + dU_att_bl.y() + dU_att_br.y();
    u3 = (j1.x()*dU_att_fl.x()+ j1.y()*dU_att_fl.y()) + (j2.x()*dU_att_fr.x()+ j2.y()*dU_att_fr.y()) + (j3.x()*dU_att_bl.x()+ j3.y()*dU_att_bl.y()) + (j4.x()*dU_att_br.x()+ j4.y()*dU_att_br.y());
    Vector_3 u(u1,u2,u3);
  /*   
    j_fl[1][3] = j1.x();    j_fl[2][3] = j1.y();
    j_fr[1][3] = j2.x();    j_fr[2][3] = j2.y();
    j_bl[1][3] = j3.x();    j_bl[2][3] = j3.y();
    j_br[1][3] = j4.x();    j_br[2][3] = j4.y();
    */
  //  Matrix u =  ~j_fl*w_fl;
    
  //  config_force(dU_att_fl, j1, u1);
  //  config_force(dU_att_fr, j2, u2);
  //  config_force(dU_att_bl, j3, u3);
  //  config_force(dU_att_br, j4, u4);
//    u_att = config_force(dU_att_fl, p.jacobian_fl()) + config_force(dU_att_fr, p.jacobian_fr()) + config_force(dU_att_bl, p.jacobian_bl()) + config_force(dU_att_br, p.jacobian_br());
//u_att = config_force(dU_att_fl, j1, u1) + config_force(dU_att_fr,  j2, u2) + config_force(dU_att_bl, j3, u3) + config_force(dU_att_br,  j4, u4);
//   u has three components here
//    u_att  = u1 + u2 + u3 +u4;
    std::cout<<"passed config space"<<std::endl;
  //  u = u_att;
    //return u_att;
 /* u.~matrix();
  j_fl.~matrix();   j_fr.~matrix();  j_bl.~matrix();  j_br.~matrix();
  w_fl.~matrix();   w_fr.~matrix();  w_bl.~matrix();  w_br.~matrix();
  p.~Listener();
*/
//u.free();
return u;
}



Point_2 att_potential_field(double dist_att, Point_2 goal, Point_2 ctrl_pt)
{
  double dist_rep;
  double att_gain = 2;
  double att_thresh = 0.5;
  double U_att;
  Point_2 dU_att;
  
  double a, b;
 
 
  a= ctrl_pt.x() - goal.x();
  b= ctrl_pt.y() - goal.y();
 Point_2 q_qgoal(a,b); 
 if (dist_att < att_thresh) 
		{
   U_att =  0.5 * att_gain * pow(dist_att,2);
//   dU_att = att_gain * q_qgoal;

   a =(double) q_qgoal.x() *att_gain;
   b =(double) q_qgoal.y() *att_gain;

   Point_2 dU_att (a,b);
 		}
  else
		{
   U_att = att_thresh * att_gain * dist_att - 0.5 * att_gain * pow(att_thresh,2);
   //dU_att = (double) (att_thresh * att_gain * q_qgoal)/dist_att;
   a = (double)(  q_qgoal.x() *att_gain*att_thresh)/dist_att;
   b =  (double)(   q_qgoal.y() *att_gain*att_thresh)/dist_att;
   Point_2 dU_att (a,b);

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

void config_force(Matrix wksp_force, Matrix jac, Matrix &cfg)
{
Matrix cfg_force;
     std::cout<<"inside config space"<<std::endl;
cfg_force = ~jac * wksp_force;

cfg = cfg_force;
//return cfg_force;

}

int main(int argc, char** argv)
{
  double dist_rep, dist_att;
  Point_2 goal(5,5);
  geometry_msgs::Twist vel;
  Listener pose;
  Vector_3 command;
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
  command=planner(pose,pgn_1,pgn_2,pgn_3,pgn_4, goal);
  std::cout<<"came out safely"<<std::endl;
  std::cout<<command<<std::endl;
  vel.linear.x = command.x();
  vel.linear.y = command.y();
  vel.linear.z = command.z();
loop_rate.sleep();
}
  
  

  return 0;
}

Listener::~Listener() {
   // Deallocate the memory that was previously reserved
   //  for this string.
  
}


