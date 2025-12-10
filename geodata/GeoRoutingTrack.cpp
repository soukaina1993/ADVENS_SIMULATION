//
// Created by cornelia.blanke on 26.06.2024.
//

#include<sstream>
#include<string>
#include <fstream>
#include "GeoRoutingTrack.h"


GeoRoutingTrack::GeoRoutingTrack()=default;

GeoRoutingTrack::GeoRoutingTrack(GeoCoord* pstart, GeoCoord* pend) : GeoTrack(pstart, pend) {}

GeoRoutingTrack::GeoRoutingTrack(GeoCoord& start, GeoCoord& end) : GeoTrack(start, end) {}

GeoRoutingTrack::~GeoRoutingTrack(){
    if (steps_created){
        for (auto step : steps)
            delete step;
    }
}

string GeoRoutingTrack::execute_curl(const string& strcurl, const short buffersize){
    // https://stackoverflow.com/questions/44610978/popen-writes-output-of-command-executed-to-cout
    char buffer[buffersize];
    string result;
    FILE* pipe = popen(strcurl.c_str(), "r");

    if (!pipe) {
        cerr << "Couldn't start command popen()." << endl;
        exit(1);
    }
    while (fgets(buffer, buffersize, pipe) != nullptr) {
        result += buffer;
    }
    pclose(pipe);
    return result;
}

void GeoRoutingTrack::get_steps() {
    compute_track();
//    cerr << "Don't use routing on demo-server!" << endl;
//    steps.resize(2);
//    steps[0] = start;
//    steps[1] = end;
}

void GeoRoutingTrack::compute_track(){
    // prepare curl call
    ostringstream curl;
    streamsize ss = curl.precision();
    curl << fixed << setprecision(10);
    //curl << "curl \"https://router.project-osrm.org/route/v1/driving/";               // demo-server (by car)
    //curl << "curl \"https://routing.openstreetmap.de/routed-car/route/v1/driving/";   // FOSSGIS demo-server (by car)
    //curl << "curl \"https://routing.openstreetmap.de/routed-foot/route/v1/driving/";  // FOSSGIS demo-server (by foot)

    // set-up of my own docker server, but use -p /opt/foot.lua
    // https://phabi.ch/2020/05/06/run-osrm-in-docker-on-windows/
    // https://gist.github.com/AlexandraKapp/e0eee2beacc93e765113aff43ec77789
    curl << "curl \"http://localhost:5000/route/v1/driving/";
    curl << start->lon << "," << start->lat << ";" << end->lon << "," << end->lat;
    curl << "?steps=false";
    curl << "&overview=simplified";     // full -> more points
    curl << "&skip_waypoints=true";
    curl << "&geometries=geojson";
    curl << "\"2>&0";                   // 2>&0 don't show error message
    curl << defaultfloat << setprecision(int(ss));
    const string strcurl = std::move(curl).str();

    string result = execute_curl(strcurl);
    //cout << result << endl;

    // process result (depends on output format)
    if (result.size() > 11 && result.substr(1, 11)==R"("code":"Ok")"){
        std::size_t found1 = result.find("coordinates");
        std::size_t found2 = result.find("]]", found1 + 14);
        string coordinates = result.substr(found1 + 14, found2 - found1 - 13);
        get_steps_from_string(coordinates);

        found1 = result.rfind("distance");
        found2 = result.find(',', found1 + 10);
        if (found2==std::string::npos)
            found2 = result.find('}', found1 + 10);
        length = stod(result.substr(found1+10, found2 - found1));
    }
    else {
        cerr << "Routing not successful." << endl;
    }
}

void GeoRoutingTrack::get_steps_from_string(const string& instring){
    std::size_t begin, mid, end, pos=0;
    steps_created = true;
    while((begin = instring.find('[', pos)) != std::string::npos) {
        mid = instring.find(',', pos);
        end = instring.find(']', pos);
        auto* step = new GeoCoord(GeoCoord::WGS84,
                                  stod(instring.substr(begin + 1, mid)), // longitude
                                  stod(instring.substr(mid + 1, end)));  // latitude
        steps.push_back(step);
        pos = end + 2;
    }
}
