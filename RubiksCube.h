#ifndef RUBIKSCUBE_H_INCLUDED
#define RUBIKSCUBE_H_INCLUDED

#include "Common.h"
#include "RubiksSide.h"
#include "RubiksColor.h"
#include "RotationDirection.h"

class RubiksCube {
private:
	int top[3][3];
	int left[3][3];
	int right[3][3];
	int front[3][3];
	int back[3][3];
	int down[3][3];

	int (*sides[6])[3][3];

	std::string result;

	void spinSide(RubiksSide side) {
		static int buffer[ 3 ];

		if (side == TOP) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = left[i][2];
			}
			for (int i = 0; i < 3; i++) {
				left[i][2] = front[0][i];
			}
			for (int i = 0; i < 3; i++) {
				front[0][i] = right[3 - i - 1][0];
			}
			for (int i = 0; i < 3; i++) {
				right[i][0] = back[2][i];
			}
			for (int i = 0; i < 3; i++) {
				back[2][3 - i - 1] = buffer[i];
			}
		} else if (side == LEFT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[i][2];
			}
			for (int i = 0; i < 3; i++) {
				down[3 - i - 1][2] = front[i][0];
			}
			for (int i = 0; i < 3; i++) {
				front[i][0] = top[i][0];
			}
			for (int i = 0; i < 3; i++) {
				top[i][0] = back[i][0];
			}
			for (int i = 0; i < 3; i++) {
				back[3 - i - 1][0] = buffer[i];
			}
		} else if (side == BACK) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[0][i];
			}
			for (int i = 0; i < 3; i++) {
				down[0][i] = left[0][i];
			}
			for (int i = 0; i < 3; i++) {
				left[0][i] = top[0][i];
			}
			for (int i = 0; i < 3; i++) {
				top[0][i] = right[0][i];
			}
			for (int i = 0; i < 3; i++) {
				right[0][i] = buffer[i];
			}
		} else if (side == RIGHT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[i][0];
			}
			for (int i = 0; i < 3; i++) {
				down[i][0] = back[3 - i - 1][2];
			}
			for (int i = 0; i < 3; i++) {
				back[i][2] = top[i][2];
			}
			for (int i = 0; i < 3; i++) {
				top[i][2] = front[i][2];
			}
			for (int i = 0; i < 3; i++) {
				front[3 - i - 1][2] = buffer[i];
			}
		} else if (side == FRONT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[2][i];
			}
			for (int i = 0; i < 3; i++) {
				down[2][i] = right[2][i];
			}
			for (int i = 0; i < 3; i++) {
				right[2][i] = top[2][i];
			}
			for (int i = 0; i < 3; i++) {
				top[2][i] = left[2][i];
			}
			for (int i = 0; i < 3; i++)
				left[2][i] = buffer[i];
		} else if (side == DOWN) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = front[2][i];
			}
			for (int i = 0; i < 3; i++) {
				front[2][i] = left[i][0];
			}
			for (int i = 0; i < 3; i++) {
				left[i][0] = back[0][3 - i - 1];
			}
			for (int i = 0; i < 3; i++) {
				back[0][i] = right[i][2];
			}
			for (int i = 0; i < 3; i++) {
				right[3 - i - 1][2] = buffer[i];
			}
		}
	}

	void spinClockwise(int side[3][3], int times, RubiksSide index) {
		static int buffer[3][3];
		static int newarray[3][3];

		if (times == 0) {
			return;
		}

		/*
		 * Transponse.
		 */
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				newarray[j][i] = side[i][j];
			}
		}
		/*
		 * Rearrange.
		 */
		for (int i = 0; i < 3; i++) {
			static int cache = 0;
			cache = newarray[i][0];
			newarray[i][0] = newarray[i][2];
			newarray[i][2] = cache;
		}

		spinSide(index);
		memcpy(buffer, newarray, sizeof(int)*3*3);

		for (int t = 1; t < times; t++) {
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < 3; i++) {
					newarray[j][i] = buffer[i][j];
				}
			}
			for (int i = 0; i < 3; i++) {
				static int cache = 0;
				cache = newarray[i][0];
				newarray[i][0] = newarray[i][2];
				newarray[i][2] = cache;
			}

			spinSide(index);

			memcpy(buffer, newarray, sizeof(int)*3*3);
		}

		memcpy(side, buffer, sizeof(int)*3*3);
	}

	double euclidean(const RubiksCube &cube) const {
		double difference = 0.0;

		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				difference += abs(top[i][j]-cube.top[i][j]);
				difference += abs(left[i][j]-cube.left[i][j]);
				difference += abs(right[i][j]-cube.right[i][j]);
				difference += abs(front[i][j]-cube.front[i][j]);
				difference += abs(back[i][j]-cube.back[i][j]);
				difference += abs(down[i][j]-cube.down[i][j]);
			}
		}

		return difference;
	}

	double colors(const RubiksCube &cube) const {
		//TODO Change array with STL maps.
		static const double coefficients[7][7] = {
			{0, 0, 0, 0, 0, 0, 0},
			{0, 1, 2, 2, 2, 2, 4},
			{0, 2, 1, 2, 4, 2, 2},
			{0, 2, 2, 1, 2, 4, 2},
			{0, 2, 4, 2, 1, 2, 2},
			{0, 2, 2, 4, 2, 1, 2},
			{0, 4, 2, 2, 2, 2, 1},
		};

		double difference = 0.0;

		/*
		 * Count matches for all sides.
		 */
		for(int s=0; s<6; s++) {
			for(int i=0; i<3; i++) {
				for(int j=0; j<3; j++) {
					/*
					 * If colors are equal calculate distance.
					 */
					difference += coefficients[(*sides[s])[1][1]][(*sides[s])[i][j]];
				}
			}
		}

		return difference;
	}

	double hausdorff(const RubiksCube &cube) const {
		long ha = 0;
		long hb = 0;
		long result = 0;

		for(int m=0; m<3; m++) {
			for(int n=0; n<3; n++) {
				int distances[] = {0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

				for(int i=0, d=0; i<3; i++) {
					for(int j=0; j<3; j++) {
						distances[d++] = abs(top[m][n]-cube.top[i][j]);
						distances[d++] = abs(left[m][n]-cube.left[i][j]);
						distances[d++] = abs(right[m][n]-cube.right[i][j]);
						distances[d++] = abs(front[m][n]-cube.front[i][j]);
						distances[d++] = abs(back[m][n]-cube.back[i][j]);
						distances[d++] = abs(down[m][n]-cube.down[i][j]);
					}
				}

				int min = distances[0];
				for(int d=0; d<54; d++) {
					if(distances[d] < min) {
						min = distances[d];
					}
				}

				if(min > ha) {
					ha = min;
				}
			}
		}

		for(int m=0; m<3; m++) {
			for(int n=0; n<3; n++) {
				int distances[] = {0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

				for(int i=0, d=0; i<3; i++) {
					for(int j=0; j<3; j++) {
						distances[d++] = abs(top[i][j]-cube.top[m][n]);
						distances[d++] = abs(left[i][j]-cube.left[m][n]);
						distances[d++] = abs(right[i][j]-cube.right[m][n]);
						distances[d++] = abs(front[i][j]-cube.front[m][n]);
						distances[d++] = abs(back[i][j]-cube.back[m][n]);
						distances[d++] = abs(down[i][j]-cube.down[m][n]);
					}
				}

				int min = distances[0];
				for(int d=0; d<54; d++) {
					if(distances[d] < min) {
						min = distances[d];
					}
				}

				if(min > hb) {
					hb = min;
				}
			}
		}

		result = std::max(ha, hb);

		return(result);
	}

	friend std::ostream& operator<< (std::ostream &out, const RubiksCube &cube);

public:
	RubiksCube() {
		reset();

		sides[0] = &top;
		sides[1] = &left;
		sides[2] = &right;
		sides[3] = &front;
		sides[4] = &back;
		sides[5] = &down;
	}

	void reset() {
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				top[i][j] = GREEN;
				left[i][j] = PURPLE;
				right[i][j] = RED;
				front[i][j] = WHITE;
				back[i][j] = YELLOW;
				down[i][j] = BLUE;
			}
		}
	}

	double compare(const RubiksCube &cube) const {
		return euclidean(cube);
	}

	void callSpin(RubiksSide side, RotationDirection direction, int numberOfTimes) {
		if (numberOfTimes < 0) {
			numberOfTimes = -numberOfTimes;
			if(direction == CLOCKWISE) {
				direction = COUNTERCLOCKWISE;
			} else if(direction == COUNTERCLOCKWISE) {
				direction = CLOCKWISE;
			}
		}

		numberOfTimes %= 4;

		if (direction == CLOCKWISE) {
			if (side == NONE) {
				/*
				* Do nothing.
				*/
			}
			if (side == TOP) {
				spinClockwise(top, numberOfTimes, TOP);
			}
			if (side == LEFT) {
				spinClockwise(left, numberOfTimes, LEFT);
			}
			if (side == RIGHT) {
				spinClockwise(right, numberOfTimes, RIGHT);
			}
			if (side == FRONT) {
				spinClockwise(front, numberOfTimes, FRONT);
			}
			if (side == BACK) {
				spinClockwise(back, numberOfTimes, BACK);
			}
			if (side == DOWN) {
				spinClockwise(down, numberOfTimes, DOWN);
			}
		}
	}

	void execute(std::string commands) {
		for(int i=0; i<commands.length(); i++) {
			callSpin((RubiksSide)commands[i], CLOCKWISE, 1);
		}
	}

	std::string shuffle(int numberOfMoves=0) {
		std::string commands = "";

		for(int i=0; i<numberOfMoves; i++) {
			switch(rand()%6) {
			case 0:
				commands+=(char)TOP;
				break;
			case 1:
				commands+=(char)LEFT;
				break;
			case 2:
				commands+=(char)RIGHT;
				break;
			case 3:
				commands+=(char)FRONT;
				break;
			case 4:
				commands+=(char)BACK;
				break;
			case 5:
				commands+=(char)DOWN;
				break;
			}
		}

		execute(commands);

		return commands;
	}

	const std::string& toString() {
		result = "";

		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(top[i][j]) + " ";
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(left[i][j]) + " ";
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(right[i][j]) + " ";
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(front[i][j]) + " ";
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(back[i][j]) + " ";
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				result += std::to_string(down[i][j]) + " ";
			}
		}

		/*
		 * Trim spaces.
		 */
		result.erase(result.size()-1, 1);
		result += '\0';

		return result;
	}

	void fromString(const char text[]) {
		std::string buffer(text);
		std::istringstream in(buffer);

		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> top[i][j];
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> left[i][j];
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> right[i][j];
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> front[i][j];
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> back[i][j];
			}
		}
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				in >> down[i][j];
			}
		}
	}
};

std::ostream& operator<< (std::ostream &out, const RubiksCube &cube) {
	for(int i=0; i<3; i++) {
		out << "      ";
		for(int j=0; j<3; j++) {
			out << cube.back[i][j] << " ";
		}
		out << std::endl;
	}

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			out << cube.left[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.top[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.right[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.down[i][j] << " ";
		}
		out << std::endl;
	}

	for(int i=0; i<3; i++) {
		out << "      ";
		for(int j=0; j<3; j++) {
			out << cube.front[i][j] << " ";
		}
		out << std::endl;
	}

	return out;
}

#endif
