//
// Created by cornelia.blanke on 25.06.2024.
//

#include "GeoCoord.h"
#include <algorithm>

GeoCoord::GeoCoord()= default;

GeoCoord::GeoCoord(const double north, const double east, const double altitude) :
            north(north), east(east), altitude_CH(altitude) {
    compute_wgs84();
}

GeoCoord::GeoCoord(GeoCoord::CRS crs, double value1, double value2, double altitude){
    if (crs == GeoCoord::CH1903plus){
        north = value1;
        east = value2;
        altitude_CH = altitude;
        compute_wgs84();
    }
    else if (crs == GeoCoord::CH1903){
        north = value1 + 2.0e6;
        east = value2 + 1.0e6;
        altitude_CH = altitude;
        compute_wgs84();
    }
    else if (crs == GeoCoord::WGS84){
        lon = value1;
        lat = value2;
        altitude_wgs84 = altitude;
        compute_CH1903plus();
    }
    else{
        cerr << "Coordinate System not implemented" << endl;
        exit(1);
    }
}

GeoCoord::GeoCoord(const GeoCoord& coord) = default;


void GeoCoord::compute_wgs84(){
    double y = (north - 26.0e5)/1e6, y2 = pow(y, 2);
    double x = (east - 12.0e5)/1e6, x2 = pow(x, 2);
    lon = (2.6779094
            + y * (4.728982 + 0.791484 * x + 0.1306 * x2 - 0.0436 * y2)) * 100.0/36.0;
    lat = (16.9023892 + 3.238272 * x - 0.270978 * y2 - 0.002528 * x2
            - x * (0.0447 * y2 + 0.0140 * x2)) * 100.0/36.0;
    altitude_wgs84 = altitude_CH + 49.55 - 12.60 * y - 22.64 * x;
}

void GeoCoord::compute_CH1903plus(){
    double lambda = 0.36 * lon - 2.67825, lambda2 = pow(lambda, 2);
    double phi =  0.36 * lat - 16.902866, phi2 = pow(phi, 2);
    north = 2600072.37 + lambda * (211455.93 - 10938.51 * phi - 0.36 * phi2 - 44.54 * lambda2);
    east = 1200147.07 + 308807.95 * phi + lambda2 * (3745.25 - 194.56 * phi)
            + phi2 * (76.63 + 119.79 * phi);
    altitude_CH = altitude_wgs84 - 49.55 + 2.73 * lambda + 6.94 * phi;
}

double GeoCoord::sqr_distance(const GeoCoord& b_coord, const bool threedim) const{
    return pow((north-b_coord.north),2) + pow((east-b_coord.east),2)
                + threedim * pow(altitude_CH - b_coord.altitude_CH, 2);
}

double GeoCoord::distance(const GeoCoord& b_coord, const bool threedim) const{
    return sqrt(sqr_distance(b_coord, threedim));
}

double distance_2d(const GeoCoord& a, const GeoCoord& b){
    return a.distance(b, false);
}

double sqr_distance_2d(const GeoCoord& a, const GeoCoord& b){
    return a.sqr_distance(b, false);
}

double distance_3d(const GeoCoord& a, const GeoCoord& b){
    return a.distance(b, true);
}

double sqr_distance_3d(const GeoCoord& a, const GeoCoord& b){
    return a.sqr_distance(b, true);
}

string GeoCoord::as_string(GeoCoord::CRS crs, const char bracket, int precision) const{
    char open, close;
    if (bracket=='('){
        open = '(';
        close = ')';
    }
    else if (bracket=='['){
        open = '[';
        close = ']';
    }
    else {
        cerr << "Option not implemented" << endl;
        exit(1);
    }
    ostringstream out;
    streamsize ss = out.precision();

    if (crs == GeoCoord::CH1903plus) {
        out << fixed << setprecision((precision >= 0) ? precision : 2 );
        out << open << north << ", " << east << close;
    }
    else if (crs == GeoCoord::CH1903) {
        out << fixed << setprecision((precision >= 0) ? precision : 2 );
        out << open << north - 2.0e6 << ", " << east - 1.0e6 << close;
    }
    else if (crs == GeoCoord::WGS84) {
        out << fixed << setprecision((precision >= 0) ? precision : 5 );
        out << open << lat << ", " << lon << close;
    }
    out << defaultfloat << setprecision(int(ss));
    return std::move(out).str();
}

void GeoCoord::display(const GeoCoord::CRS crs) const{
    streamsize ss = cout.precision();
    if (crs == GeoCoord::CH1903plus) {
        cout << fixed << setprecision(2);
        cout << "North " << north << "\tEast " << east << "\t\tAltitude " << altitude_CH << "\n";
    }
    else if (crs == GeoCoord::CH1903) {
        cout << fixed << setprecision(2);
        cout << "North " << north - 2.0e6 << "\tEast " << east - 1.0e6 << "\t\tAltitude " << altitude_CH << "\n";
    }
    else if (crs == GeoCoord::WGS84) {
        cout << fixed << setprecision(6);
        cout << "Latitude " << lat << "\tLongitude " << lon << "\tAltitude " << altitude_wgs84 << "\n";
    }
    else {
        cerr << "Coordinate System not implemented" << endl;
        exit(1);
    }
    cout << defaultfloat << setprecision(int(ss)) << endl;
}

void GeoCoord::display() const{
    streamsize ss = cout.precision();
    cout << fixed << setprecision(2);
    cout << "Coordinate System CH1903+\n";
    cout << "North " << north << "\tEast " << east << "\t\tAltitude " << altitude_CH << "\n";
    cout << fixed << setprecision(6);
    cout << "Coordinate System WGS84\n";
    cout << "Latitude " << lat << "\tLongitude " << lon << "\tAltitude " << altitude_wgs84 << "\n";
    cout << defaultfloat << setprecision(int(ss)) << endl;
}
