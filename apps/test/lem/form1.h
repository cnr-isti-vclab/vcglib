/****************************************************************************
** Form interface generated from reading ui file 'form1.ui'
**
** Created: Mon Aug 2 16:35:19 2004
**      by: The User Interface Compiler ($Id: form1.h,v 1.1 2004-08-04 20:59:13 pietroni Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef FORM1_H
#define FORM1_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;

class Form1 : public QDialog
{
    Q_OBJECT

public:
    Form1( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~Form1();


protected:

protected slots:
    virtual void languageChange();

};

#endif // FORM1_H
