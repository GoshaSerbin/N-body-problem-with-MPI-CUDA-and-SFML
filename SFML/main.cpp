#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <SFML/Graphics.hpp>


using namespace sf;

typedef float T;

#define  MAX(a,b) ((a)>(b)) ? (a) : (b)

const T PI = 3.14159265358979;

const T k = 0.9;

const T G = 40;
const T eps = 20;

std::random_device rd;
std::mt19937 gen(rd());

#define MAX_VAL 1e13
std::uniform_real_distribution<> distribution(0, MAX_VAL);

T fRand(T fMin, T fMax)
{
	T f = distribution(gen) / MAX_VAL;
	return fMin + f * (fMax - fMin);
}

struct point
{

	T coord[2];
	T vel[2];
	T mass;


	point& operator= (const point& p)
	{
		mass = p.mass;
		for (int i = 0; i < 2; ++i) {
			coord[i] = p.coord[i];
			vel[i] = p.vel[i];
		}
		return *this;
	}

	point() {
		coord[0] = 0; coord[1] = 0;
		vel[0] = 0; vel[1] = 0;
		mass = 1;
	}
	point(T COORD[], T VEL[], T MASS) {
		coord[0] = COORD[0]; coord[1] = COORD[1];
		vel[0] = VEL[0]; vel[1] = VEL[1];
		mass = MASS;
	}

	T getRadius() {
		return 10 * sqrt(mass);
	}

};

std::ostream& operator<< (std::ostream& out, const point& point)
{
	out << point.coord[0] << " " << point.coord[1];
	return out;
}

T* f(point* points, int N, T* accel, T* result) {
	//y = {rx,ry,vx,vy}
	for (int i = 0; i < N * 4; i += 4) {
		result[i] = points[i / 4].vel[0];
		result[i + 1] = points[i / 4].vel[1];
		result[i + 2] = accel[i / 2];
		result[i + 3] = accel[i / 2 + 1];
	}

	return result;
}

T dnorm(const T* vec1, const T* vec2)
{
	T sum = 0;
	for (int i = 0; i < 2; ++i)
	{
		sum += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
	}
	return sqrt(sum);
}

T norm(T vec[]) {
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
}
void calculateAccelerations(const int N, const point* points, T* accel)
{
	for (int pi = 0; pi < N; ++pi) {
		point p = points[pi];
		for (int dim = 0; dim < 2; ++dim)
		{
			T a = 0;
			for (int i = 0; i < N; ++i)
			{
				T denominator = pow(dnorm(p.coord, points[i].coord), 3);
				a += points[i].mass * (p.coord[dim] - points[i].coord[dim]) / (MAX(denominator, eps));
			}

			accel[2 * pi + dim] = -G * a;
		}
	}

}



void addPoint(point pnt, point* points, int& real_number) {
	points[real_number] = pnt;
	real_number++;
}

void draw(point* points, RenderWindow& window, int real_number) {//can be updated
	
	
	for (int i = 0; i < real_number; ++i) {
		CircleShape point_shape(10* sqrt(points[i].mass), 20);
		T radius = point_shape.getRadius();
		point_shape.setFillColor(Color::Green);
		point_shape.setPosition(points[i].coord[0] - radius, points[i].coord[1] - radius);
		window.draw(point_shape);
	}
}

int main()
{
    const int SCREEN_WIDTH = 960;
    const int SCREEN_HEIGHT = 750;
    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "PROJECT_SFML"/*,Style::Fullscreen*/);


	int N = 500;
	int real_number = 0;
	Clock clock;
	T time = -0.01;
	T tau;
	
	point* points = new point[N];
	bool* isOnCollision = new bool[N];

	point* ps1 = new point[N];

	T* k1 = new T[N * 4];
	T* k2 = new T[N * 4];

	T* accel = new T[N * 2];
	
	T coord[2] = { SCREEN_WIDTH/2, SCREEN_HEIGHT/2 };
	T vel[2] = { 0,0 };
	T mass = 100;
	point pnt(coord, vel, mass);
	addPoint(pnt, points, real_number);

	bool isCollising = false;

	while (window.isOpen()) {

		Vector2i mouse_pos = Mouse::getPosition(window);
		//int rect_width = window.getSize().x;
		//int rect_height = window.getSize().y;

		Event event;

		while (window.pollEvent(event)) {

			if (event.type == Event::Closed)	window.close();

			if (event.type == Event::MouseButtonPressed && event.key.code == Mouse::Left) {
				//distance = hypot(points[0].coord[0] - points[1].coord[0], points[0].coord[1] - points[1].coord[1]);
				//if (distance < st_shape.getRadius()) {
				//isCreatingParticle = true;
				
				//}

			}

			if (event.type == Event::MouseButtonReleased && event.key.code == Mouse::Left) {
				
				T mass = fRand(0.5, 1);
				T coord[2] = { mouse_pos.x, mouse_pos.y};
				T vel[2] = { fRand(-1,1), fRand(-1,1) };
				point pnt(coord, vel, mass);
				addPoint(pnt, points, real_number);
				std::cout << "position = " << points[real_number - 1] << "\n";
				std::cout << "real_number ="<<real_number<<" \n";
				
			}



		}
		//std::cout << points[0] << "\n";
		//if (isSelected) {
		//	//st_shape.setPosition(Vector2f(mouse_pos.x - st_shape.getRadius(), mouse_pos.y - st_shape.getRadius()));
		//}
		tau = 50*(clock.getElapsedTime().asSeconds() - time);
		time = clock.getElapsedTime().asSeconds();

		//movement
		{

			const T tau05 = tau / 2.0;

			calculateAccelerations(real_number, points, accel);


			f(points, N, accel, k1);
			for (int pi = 0; pi < real_number; ++pi)
			{
				ps1[pi] = points[pi];


				for (int j = 0; j < 2; ++j) {
					points[pi].coord[j] = ps1[pi].coord[j] + tau05 * k1[4 * (pi)+j];
					points[pi].vel[j] = ps1[pi].vel[j] + tau05 * k1[4 * (pi)+2 + j];
				}

			}

			calculateAccelerations(real_number, points, accel);

			f(points, real_number, accel, k2);


			for (int pi = 0; pi < real_number; ++pi) {

				for (int j = 0; j < 2; ++j) {
					points[pi].coord[j] = ps1[pi].coord[j] + tau * k2[4 * (pi)+j];
					points[pi].vel[j] = ps1[pi].vel[j] + tau * k2[4 * (pi)+2 + j];
				}

			}

		}



		//for (int i = 0; i < real_number; ++i) {
		//	points[i].coord[0] += points[i].vel[0] * tau;
		//	points[i].coord[1] += points[i].vel[1] * tau;
		//}

				//collision
		{
			for (int i = 0; i < real_number; ++i) {
				point pnt_i = points[i];
				for (int j = i + 1; j < real_number; ++j) {
					point pnt_j = points[j];
					if (dnorm(pnt_i.coord, pnt_j.coord) < pnt_i.getRadius() + pnt_j.getRadius()) {
						std::cout << "Collision!\n";
						
							//points[i].coord[0] -= points[i].vel[0] * 2 *tau;
							//points[i].coord[1] -= points[i].vel[1] * 2 * tau;
							//points[j].coord[0] -= points[j].vel[0] * 2 * tau;
							//points[j].coord[1] -= points[j].vel[1] * 2 * tau;
						//if (!isOnCollision[i] & !isOnCollision[j]) {
						//	isOnCollision[i] = true; isOnCollision[j] = true;
							T teta_i = atan2f(pnt_i.vel[1], pnt_i.vel[0]);
							T teta_j = atan2f(pnt_j.vel[1], pnt_j.vel[0]);
							T phi = atan2f(pnt_j.coord[1] - pnt_i.coord[1], pnt_j.coord[0] - pnt_i.coord[0]);
							T m_i = pnt_i.mass;
							T m_j = pnt_j.mass;
							T v_i = norm(pnt_i.vel);
							T v_j = norm(pnt_j.vel);
							points[i].vel[0] = (v_i * cos(teta_i - phi) * (m_i - m_j) + 2 * m_j * v_j * cos(teta_j - phi)) * cos(phi) / (m_i + m_j)
								+ v_i * sin(teta_i - phi) * cos(phi + PI / 2);
							points[i].vel[1] = (v_i * cos(teta_i - phi) * (m_i - m_j) + 2 * m_j * v_j * cos(teta_j - phi)) * sin(phi) / (m_i + m_j)
								+ v_i * sin(teta_i - phi) * sin(phi + PI / 2);
							points[j].vel[0] = (v_j * cos(teta_j - phi) * (m_j - m_i) + 2 * m_i * v_i * cos(teta_i - phi)) * cos(phi) / (m_j + m_i)
								+ v_j * sin(teta_j - phi) * cos(phi + PI / 2);
							points[j].vel[1] = (v_j * cos(teta_j - phi) * (m_j - m_i) + 2 * m_i * v_i * cos(teta_i - phi)) * sin(phi) / (m_j + m_i)
								+ v_j * sin(teta_j - phi) * sin(phi + PI / 2);
						//}
							points[i].vel[0] *= sqrt(k);
							points[i].vel[1] *= sqrt(k);
							points[j].vel[0] *= sqrt(k);
							points[j].vel[1] *= sqrt(k);

					}
					//else {
					//	isOnCollision[i] = false; isOnCollision[j] = false;
					//}
				}
			}
		}

		for (int i = 0; i < real_number; ++i) {
			if (points[i].coord[0]+ points[i].getRadius() > SCREEN_WIDTH || points[i].coord[0] - points[i].getRadius() < 0) {
				points[i].vel[0] *= -1;
			}
			if (points[i].coord[1] + points[i].getRadius() > SCREEN_HEIGHT || points[i].coord[1] - points[i].getRadius() < 0) {
				points[i].vel[1] *= -1;
			}
		}


		window.clear();

		draw(points, window, real_number);

		window.display();
		//sleep(seconds(1));
	}


	delete[] k1;
	delete[] k2;
	delete[] accel;
	delete[] ps1;

	return 0;

}