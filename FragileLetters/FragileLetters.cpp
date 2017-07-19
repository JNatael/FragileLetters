//Defines the entry point for the console application.
#include <iostream>
#include <vector>
#include <stack>
#include <math.h>
#include <algorithm>
using namespace std;

/*
This function reads a line of int values into a vector function and returns that
vector.
*/
vector<int> readlineofints(int count) {
    // Input values
    vector<int> linevalues(count);
    for (int j = 0; j < count; j++) {
        int val;
        cin >> val;
        linevalues[j] = val;
    }
    return linevalues; // Return line values as a vector
}

vector<double> readlineofdoubles(int count) {
    // Input values
    vector<double> linevalues(count);
    for (int j = 0; j < count; j++) {
        double val;
        cin >> val;
        linevalues[j] = val;
    }
    return linevalues; // Return line values as a vector
}

struct Point {
    double x;
    double y;

    //Add back default constructor
    Point() = default;   

    //Constructor
    Point(double x_in, double y_in) :x(x_in), y(y_in) {};

    //Copy constructor
    Point(const Point& input) :x(input.x), y(input.y) {};
};

//Pointer comparison operator
bool operator==(const Point &a, const Point &b) {
    return(a.x == b.x && a.y == b.y);
}


//Taken from http://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/
// A globle point needed for  sorting points with reference
// to  the first point Used in compare function of qsort()
Point p0;

// A utility function to find next to top in a stack
Point nextToTop(stack<Point> &S) {
    Point p = S.top();
    S.pop();
    Point res = S.top();
    S.push(p);
    return res;
}

// A utility function to swap two points
void swap(Point &p1, Point &p2) {
    Point temp = p1;
    p1 = p2;
    p2 = temp;
}

// A utility function to return square of distance
// between p1 and p2
double distSq(Point p1, Point p2) {
    return (p1.x - p2.x)*(p1.x - p2.x) +
        (p1.y - p2.y)*(p1.y - p2.y);
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point p, Point q, Point r) {
    long val = (q.y - p.y) * (r.x - q.x) -
        (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // colinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void *vp1, const void *vp2) {
    Point *p1 = (Point *)vp1;
    Point *p2 = (Point *)vp2;

    // Find orientation
    int o = orientation(p0, *p1, *p2);
    if (o == 0)
        return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

    return (o == 2) ? -1 : 1;
}

// Finds convex hull of a set of n points.
vector<Point> convexHull(Point points[], int n) {
    // Find the bottommost point
    double ymin = points[0].y;
    long min = 0;
    for (int i = 1; i < n; i++) {
        double y = points[i].y;

        // Pick the bottom-most or chose the left
        // most point in case of tie
        if ((y < ymin) || (ymin == y &&
            points[i].x < points[min].x))
            ymin = points[i].y, min = i;
    }

    // Place the bottom-most point at first position
    swap(points[0], points[min]);

    // Sort n-1 points with respect to the first point.
    // A point p1 comes before p2 in sorted output if p2
    // has larger polar angle (in counterclockwise
    // direction) than p1
    p0 = points[0];
    qsort(&points[1], n - 1, sizeof(Point), compare);

    // If two or more points make same angle with p0,
    // Remove all but the one that is farthest from p0
    // Remember that, in above sorting, our criteria was
    // to keep the farthest point at the end when more than
    // one points have same angle.
    int m = 1; // Initialize size of modified array
    for (int i = 1; i<n; i++) {
        // Keep removing i while angle of i and i+1 is same
        // with respect to p0
        while (i < n - 1 && orientation(p0, points[i],
            points[i + 1]) == 0)
            i++;


        points[m] = points[i];
        m++;  // Update size of modified array
    }


    vector<Point> output_points;
    // If modified array of points has less than 3 points,
    // convex hull is not possible; return empty array
    if (m < 3) return output_points;

    // Create an empty stack and push first three points
    // to it.
    stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);

    // Process remaining n-3 points
    for (int i = 3; i < m; i++) {
        // Keep removing top while the angle formed by
        // points next-to-top, top, and points[i] makes
        // a non-left turn
        while (orientation(nextToTop(S), S.top(), points[i]) != 2)
            S.pop();
        S.push(points[i]);
    }


    // Now stack has the output points, print contents of stack
    while (!S.empty()) {
        Point p = S.top();
        output_points.push_back(p);
        S.pop();
    }
    return output_points;
}


Point vector_add(Point &a, Point &b) {
    return(Point((a.x + b.x), (a.y + b.y)));
}

Point vector_subtract(Point &a, Point &b) {
    return(Point((a.x - b.x), (a.y - b.y)));
}

double vector_dot_product(Point &a, Point &b) {
    return(a.x*b.x + a.y*b.y);
}

double vector_length(Point &a) {
    Point vector_btw = a;
    return(sqrt(vector_btw.x*vector_btw.x + vector_btw.y*vector_btw.y));
}

double vector_length(Point &start, Point &end) {
    Point vector_btw = vector_subtract(end, start);
    return(vector_length(vector_btw));
}

Point vector_mult(Point &a, double mult) {
    return(Point(a.x*mult, a.y*mult));
}

int find_in_orig_points(const Point (&points)[1000],const Point &input_point, const int &n) {
    //Find this point on the original polygon
    int k;
    for (k = 0; k < n; k++) {
        if (points[k] == input_point) {
            return(k);
        }
    }
    return -1;
}

/*Note: This assumes the vectors have the same start point
*/
Point vector_proj(Point &vector_start, Point &vector_end, Point &proj_pt) {
    Point main_vec = vector_subtract(vector_end, vector_start);
    Point proj_vec = vector_subtract(proj_pt, vector_start);
    double scaling = vector_dot_product(proj_vec, main_vec) / vector_dot_product(main_vec, main_vec);
    Point result = vector_mult(main_vec, scaling);
    return(vector_add(result,vector_start));
}

/*This function checks to see if the input point is on the input line segment
It assumes the point is on the line at least already (if not in the segment)*/
bool point_on_line_in_line_segment(const Point &line_point_1,const Point &line_point_2, const  Point &point_to_check) {
    bool x_in_segment = false, y_in_segment = false;
    //If x is between the line segment points (check in both directions, as we don't know whether start or end point is bigger)
    if (point_to_check.x <= line_point_2.x && point_to_check.x >= line_point_1.x) { x_in_segment = true; }
    if (point_to_check.x <= line_point_1.x && point_to_check.x >= line_point_2.x) { x_in_segment = true; }

    //Same check for y
    if (point_to_check.y <= line_point_2.y && point_to_check.y >= line_point_1.y) { y_in_segment = true; }
    if (point_to_check.y <= line_point_1.y && point_to_check.y >= line_point_2.y) { y_in_segment = true; }

    //IF both x and y were in the segment return true, otherwise return false
    if (x_in_segment && y_in_segment) {
        return true;
    }
    else {
        return false;
    }
}


int main() {
    std::ios_base::sync_with_stdio(false);

    // get test case count
    int t;
    std::cin >> t;

    //! loop over all the test cases
    for (int i = 1; i <= t; i++) {
        // Read in params
        vector<int> params = readlineofints(1);
        long n = params[0];

        //Loop over inputs
        Point points[1000], hull_point_set[1000];
        for (int j = 0; j < n; j++) {
            vector<double> this_point = readlineofdoubles(2);
            Point x = { this_point[0],this_point[1] };
            points[j] = x;
            hull_point_set[j] = x;
        }        

        //Code taken from similar code at https://gist.github.com/listochkin/1200393
        double area = 0.0, c_x = 0.0, c_y = 0.0;
        for (size_t j = 0; j < n; j++) {
            size_t k = (j + 1) % n;
            double det = points[j].x * points[k].y - points[j].y * points[k].x;
            area += det;
            c_x += (points[j].x + points[k].x)*det;
            c_y += (points[j].y + points[k].y)*det;
        }
        area = fabs(area)/2; //Area of the hull
        Point center_of_mass(fabs(c_x) / (6*area), fabs(c_y) / (6*area));

        vector<Point> hull_points;
        //Get convex hull points
        hull_points = convexHull(hull_point_set, n);

        //Iterate over convex hull
        int possible_orientations = 0;
        for (size_t j = 0; j < hull_points.size(); j++) {

            //If the convex hull edge matches one of the original edges it's a possible solution
            int orig_point_index = find_in_orig_points(points, hull_points[j], n);
            if (orig_point_index != -1) {
                Point next_hull_point = hull_points[(j + 1) % hull_points.size()];
                int next_edge_point_index = find_in_orig_points(points, next_hull_point, n);
                bool next_edge_point_is_next = ((next_edge_point_index == (orig_point_index + 1) % (n)));
                int previous_edge = orig_point_index - 1;
                if (previous_edge < 0) { previous_edge = n-1; }
                if (next_edge_point_is_next || next_edge_point_index == previous_edge) { //If the next hull point is an adjacent orig point
                    //If it's a possible solution, verify the center would be above it when it's down
                    Point projected_point = vector_proj(hull_points[j], next_hull_point, center_of_mass);
                    if(point_on_line_in_line_segment(hull_points[j],next_hull_point,projected_point)){
                        possible_orientations++; //If it is, increment the orientations possible
                    }
                }
            }
        }

        //Print output
        std::cout << "Case #" << i << ": " << possible_orientations << std::endl;

    }
    return 0;
}