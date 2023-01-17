#include <bits/stdc++.h>

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 10000;
const double small_radius = 5000;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 1000;
constexpr double PI = 3.14159265358979323846264338;
const int num_channel = 6;
const int num_power = 4;
const int num_BS = 100;
const double pass_loss_exponent = 2.5;
const double POWER = 1.0;
const double theta = pow(10, 0.4);

//string filename = "NOMA_4_BS2.txt";
//出力ファイルを開く
ofstream outputfile;

struct terminal {
    pair<double, double> pos;
    //double to_BS[num_BS];
    double nearest_dst;
    int nearest_BS;
    //double gain;
    double coef[num_BS];
    bool state; //True: in , False: out
};

struct BS {
    int channel[num_channel];
    int power[num_power];
};

//int num_terminal = 10000;
//vector<terminal> T_vec(num_terminal);


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//擬似乱数発生
double my_rand(double Min, double Max) {
    mt19937 mt{ std::random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    double g = distr(mt);
    return g;
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    double r = radius * sqrt(my_rand(0, 1));
    double theta = my_rand(-PI, PI);
    double x = r * cos(theta);
    double y = r * sin(theta);
    pair<double, double> pos = make_pair(x, y);
    return pos;
}

//2点間の距離を計算
double cal_dst(pair<double, double> pos1, pair<double, double> pos2) {
    return sqrt((pos1.first - pos2.first)*(pos1.first - pos2.first) + (pos1.second - pos2.second)*(pos1.second - pos2.second));
}


//正規分布乱数発生
double gauss_rand(double mu, double sig) {
    //mt19937 mt{ std::random_device{}() };
    random_device rd;
    default_random_engine eng(rd());
    normal_distribution<> dist(mu, sig);
    //uniform_real_distribution<double> distr(Min, Max);
    double g = dist(eng);
    return g;
}

//初期化関数
void initialization(vector<terminal>& T_vec) {
    //BSの座標を設定
    vector<pair<double, double>> BS_pos(num_BS);
    for (int i = 0; i < num_BS; i++) BS_pos.at(i) = coordinate();
    pair<double, double> origin = make_pair(0, 0); //原点座標
    
    //端末情報を初期化
    for (int i = 0; i < T_vec.size(); i++) {
        T_vec.at(i).pos = coordinate();
        //内円の内部にあればTrue
        if (cal_dst(origin, T_vec.at(i).pos) < small_radius) T_vec.at(i).state = true;
        else T_vec.at(i).state = false;
        
        T_vec.at(i).nearest_dst = inf;
        for (int j = 0; j < num_BS; j++) {
            //端末から各BSまでの距離を保存
            double dst = cal_dst(T_vec.at(i).pos, BS_pos.at(j));
            //端末から最も近いBSまでの距離を保存
            if (dst < T_vec.at(i).nearest_dst) {
                T_vec.at(i).nearest_dst = dst;
                T_vec.at(i).nearest_BS = j;
            }
            //端末-基地局間のフェージング係数を設定
            double H = gauss_rand(0, 1);
            T_vec.at(i).coef[j] = H / pow(dst, pass_loss_exponent);
        }
        
    }
    cout << "Initialization complete" << endl;
}


void simulation (vector<terminal>& T_vec, double Pr_ac) {
    if (T_vec.size() == 0) return;
    double thp_sum = 0;
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int t  = 0; t < end_time; t++) {
        vector<double> SI(num_BS, 0);
        //アクセスする端末を決定
        for (int i = 0; i < T_vec.size(); i++) {
            double r = my_rand(0, 1);
            if (r < Pr_ac) {
                for (int j = 0; j < num_BS; j++) {
                    SI.at(j) += T_vec.at(i).coef[j];
                }
            } else continue;
        }
        
        for (int i = 0; i < T_vec.size(); i++) { //修正
            if (!T_vec.at(i).state) continue;
            int x = T_vec.at(i).nearest_BS;
            double SINR = pow(T_vec.at(i).nearest_dst, 2 * pass_loss_exponent) * T_vec.at(i).coef[x] / SI.at(x);
            if (SINR > theta) thp_sum++;
        }
                
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    double throughput = thp_sum / end_time;
    cout << throughput << endl << endl;
    outputfile << throughput;
}



int main() {
    string filename = "SPALOHA.txt";
    outputfile.open(filename);
    int num_terminal = 0;
    for (int i = 0; i < 10000; i++) {
        vector<terminal> T_vec(num_terminal);
        initialization(T_vec);
        cout << num_terminal << " terminals" << endl;
        outputfile << num_terminal;
        for (int j = 0; j < 3; j++) {
            outputfile << " ";
            double Pr_ac = 0.001 * pow(10, j);
            simulation(T_vec, Pr_ac);
        }
        outputfile << endl;
        cout << endl << endl;
        num_terminal += 2000;
    }
    outputfile.close();
}
