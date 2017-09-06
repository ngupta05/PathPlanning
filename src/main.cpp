#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <stdlib.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

const double MPH_TO_MPS = 1609.34 / 3600;
double ref_v = 49 * MPH_TO_MPS; // reference vel in m/s
bool lane_change_engaged = false;
int target_ego_lane = -1;
double target_ego_v = ref_v;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            int prev_size = previous_path_x.size();

            const double PRED_STEPS_HORIZON = 50;
            double car_v = car_speed * MPH_TO_MPS;

          	json msgJson;

            vector<double> ptsx;
            vector<double> ptsy; 

            double ref_x = car_x;
            double ref_y = car_y;
            int cur_car_lane = car_d / 4;
            int end_car_lane = end_path_d / 4;
            double end_car_s = end_path_s;
            double end_car_v = car_v;

            double ref_yaw = deg2rad(car_yaw);

            if (prev_size < 2) {
              double speed = std::max(car_v, 1.0);
              double prev_car_x = car_x - speed * 0.02 * cos(car_yaw);
              double prev_car_y = car_y - speed * 0.02 * sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);

              end_car_lane = car_d / 4;
              end_car_s = car_s;
            } else {
              ref_x = previous_path_x[prev_size - 1];
              ref_y = previous_path_y[prev_size - 1];
              double prev_ref_x = previous_path_x[prev_size - 2];
              double prev_ref_y = previous_path_y[prev_size - 2];
              ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);
              ptsx.push_back(prev_ref_x);
              ptsx.push_back(ref_x);
              ptsy.push_back(prev_ref_y);
              ptsy.push_back(ref_y);

              end_car_v = sqrt(std::pow(ref_x - prev_ref_x, 2) + std::pow(ref_y - prev_ref_y, 2)) / 0.02;
            }
            
            if (lane_change_engaged && 
                cur_car_lane == end_car_lane &&
                cur_car_lane == target_ego_lane) {
              // lane change completed
              lane_change_engaged = false;
            }

            
            if (!lane_change_engaged) {
              vector<double> leading_speed = {ref_v, ref_v, ref_v};
              vector<double> leading_dist = {100000, 100000, 100000};
              vector<bool> col = {false, false, false};
              vector<bool> pot_col = {false, false, false};

            
              // sensor fusion data: id, x, y, vx, vy, s, d
              for (int i = 0; i < sensor_fusion.size(); i++) {
                int sd = sensor_fusion[i][6];
                double ss = sensor_fusion[i][5];
                double svx = sensor_fusion[i][3];
                double svy = sensor_fusion[i][4];
                double sv = sqrt(svx * svx + svy * svy) * MPH_TO_MPS;
                int slane = sd / 4;
                if (ss >= car_s && ss < (end_car_s + 20) && ss < leading_dist[slane]){
                  leading_dist[slane] = ss;
                  leading_speed[slane] = std::min(ref_v, sv);
                }

                if (slane != end_car_lane && ss > car_s - 20 && ss < end_car_s + 20) {
                  pot_col[slane] = true;
                }

                // for collision, following holds:
                // end position of ego and other car would be same at some time t
                // from end_path
                // assuming constant velocity for both
                // ss + (prev_size * 0.02 + t) * sv = end_car_s + end_car_v * t;
                double tc = (ss + prev_size * 0.02 * sv - end_car_s) /
                            std::max(end_car_v - sv, 0.01);
                bool too_close = (ss > car_s - 10 && ss < end_car_s + 10);

                if (slane != end_car_lane && (too_close || (tc > -1 && tc <= 2))) {
                  col[slane] = true;
                }
              }

              vector<int> next_lanes = {end_car_lane};

              if (end_car_lane != 0)
                next_lanes.push_back(end_car_lane - 1);
              if (end_car_lane != 2)
                next_lanes.push_back(end_car_lane + 1);

              int min_cost = 1000000;
            
              for (int i = 0; i < next_lanes.size(); i++) {
                double cost = 0;
                int l = next_lanes[i];
                double speed = leading_speed[l];

                // speed cost
                cost += (ref_v - speed) * (ref_v - speed);

                // collision cost
                if (col[l]) {
                  cost += 10000;
                }

                if (pot_col[l]) {
                  cost += 1000;
                }

                if (l == end_car_lane) {
                  bool prefer_lane_change = (rand() % 100 == 0);
                  if (prefer_lane_change) {
                    cost += 5; // tie break cases when leading v is same
                  }
                }

                if (cost < min_cost) {
                  min_cost = cost;
                  target_ego_v = speed;
                  target_ego_lane = l;
                }
              }
            }
            
            if (target_ego_lane != cur_car_lane) {
              lane_change_engaged = true;
            }
           
            double next_ego_v = end_car_v + (target_ego_v - end_car_v)/25;
            
            vector<vector<double>> next_mp;
            double wp_delta = 30;
            for (int i = 0; i < 3; i++) {
              double delta = (i + 1) * wp_delta;
              vector<double> wp = getXY(end_car_s + delta, 2 + 4 * target_ego_lane,
                map_waypoints_s, map_waypoints_x, map_waypoints_y);
              next_mp.push_back(wp);
              ptsx.push_back(wp[0]);
              ptsy.push_back(wp[1]);
            }

            // convert to car local frame ref
            for (int i = 0; i < ptsx.size(); i++) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
              ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
            }
            
            tk::spline s;
            s.set_points(ptsx, ptsy);

            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt(target_x * target_x + target_y * target_y);

          	vector<double> next_x_vals = previous_path_x;
          	vector<double> next_y_vals = previous_path_y;

            for (int i = 1; i <= PRED_STEPS_HORIZON - previous_path_x.size();
                i++) {
              double N = (target_dist/(0.02 * next_ego_v));
              double x_point = i * target_x / N;
              double y_point = s(x_point);

              double x_ref = x_point;
              double y_ref = y_point;

              x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
              y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw); 

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }
            
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































