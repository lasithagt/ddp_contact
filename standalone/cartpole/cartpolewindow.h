#ifndef CARTPOLEWINDOW_H
#define CARTPOLEWINDOW_H

#include <QMainWindow>
#include <vector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    using Scalar = double;
    explicit MainWindow(Scalar anglef, Scalar positionf, Scalar dt, QWidget *parent = 0);
    ~MainWindow();

public slots:
    void update(const QVector<Scalar> &x, const QVector<Scalar> &u);

private:
    Scalar anglef_;
    Scalar positionf_;
    Scalar dt_;
    Ui::MainWindow *ui;
};

#endif // CARTPOLEWINDOW_H
