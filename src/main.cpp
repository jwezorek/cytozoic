#include "main_window.hpp"
#include <QApplication>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <qdebug.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = CGAL::Point_2<Kernel>;

int main(int argc, char *argv[]) {


    Point_2 p(1.0, 2.0), q(4.0, 6.0);

    qDebug() << "Point p: " << p.x() << "," << p.y();
    qDebug() << "Point q: " << q.x() << "," << q.y();
    qDebug() << "Distance between p and q: " << CGAL::sqrt(CGAL::squared_distance(p, q));

    QApplication app(argc, argv);
    
    cz::main_window window;
    window.show();

    return app.exec();
}