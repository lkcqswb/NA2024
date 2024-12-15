#include "BSpline.hpp"
#include "ppForm.hpp"

using json = nlohmann::json;
using namespace std;
int main() {
    json j;
    std::ifstream file("template.json"); 
    if (!file.is_open()) {
            std::cerr << "无法打开文件" << std::endl;
            return 1;
    }
    j = json::parse(file);
    file.close();
    BSpline b(j);
    ppForm p1(j,1);
    ppForm p3(j,3);
    ppForm p2(j,2);
    double e=1e-6;
    double t=4;
    cout << (b.get_value(t-e)[0] + b.get_value(t+e)[0] - 2*b.get_value(t)[0]) / (e*e) << endl;
    cout << (p1.get_value(t-e)[0] + p1.get_value(t+e)[0] - 2*p1.get_value(t)[0]) / (e*e) << endl;
    cout << (p2.get_value(t-e)[0] + p2.get_value(t+e)[0] - 2*p2.get_value(t)[0]) / (e*e) << endl;
    cout << (p3.get_value(t-e)[0] + p3.get_value(t+e)[0] - 2*p3.get_value(t)[0]) / (e*e) << endl;

}