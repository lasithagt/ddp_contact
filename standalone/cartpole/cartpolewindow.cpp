#include "cartpolewindow.h"
#include "ui_cartpolewindow.h"

MainWindow::MainWindow(Scalar anglef, Scalar positionf, Scalar dt, QWidget *parent) :
    QMainWindow(parent),
    anglef_(anglef), positionf_(positionf), dt_(dt),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Angle
    ui->angle_plot->addGraph();
    ui->angle_plot->addGraph();
    ui->angle_plot->legend->setVisible(true);
    ui->angle_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->angle_plot->graph(0)->setName("Target");
    ui->angle_plot->graph(1)->setName("Trajectory");
    ui->angle_plot->xAxis->setLabel("Time (s)");
    ui->angle_plot->yAxis->setLabel("Angle (rad)");
    ui->angle_plot->graph(0)->setPen(QPen(Qt::red));
    ui->angle_plot->rescaleAxes(true);
    ui->angle_plot->yAxis->setRange(-2.0 * M_PI, 2.0 * M_PI);

    // Position
    ui->position_plot->addGraph();
    ui->position_plot->addGraph();
    ui->position_plot->legend->setVisible(true);
    ui->position_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->position_plot->graph(0)->setName("Target");
    ui->position_plot->graph(1)->setName("Trajectory");
    ui->position_plot->xAxis->setLabel("Time (s)");
    ui->position_plot->yAxis->setLabel("Position (m)");
    ui->position_plot->graph(0)->setPen(QPen(Qt::red));
    ui->position_plot->rescaleAxes(true);
    ui->position_plot->yAxis->setRange(-5.0, 5.0);

    // Control
    ui->control_plot->addGraph();
    ui->control_plot->graph(0)->setName("Control");
    ui->control_plot->xAxis->setLabel("Time (s)");
    ui->control_plot->yAxis->setLabel("Control");
    ui->control_plot->rescaleAxes(true);
}

void MainWindow::update(const QVector<Scalar> &x, const QVector<Scalar> &u)
{
    static Scalar time = 0.0;
    time += dt_;
    ui->position_plot->graph(0)->addData(time, positionf_);
    ui->position_plot->graph(1)->addData(time, x[0]);
    ui->angle_plot->graph(0)->addData(time, anglef_);
    ui->angle_plot->graph(1)->addData(time, x[2]);
    ui->control_plot->graph(0)->addData(time, u[0]);
    ui->position_plot->rescaleAxes(true);
    ui->position_plot->replot();
    ui->angle_plot->rescaleAxes(true);
    ui->angle_plot->replot();
    ui->control_plot->rescaleAxes(true);
    ui->control_plot->replot();
}

MainWindow::~MainWindow()
{
    delete ui;
}
