#ifndef GETOPT_H
#define GETOPT_H

#include <QString>
#include <QStringList>
#include <QMap>

/* Example usage: 

  QString file1, file2;
  QString o_gamma = "10";
  QString o_scale = "0.7";

  GetOpt opt(argc, argv);
  opt.addArgument("img1", "first image", &file1);
  opt.addArgument("img2", "second image", &file2);
  opt.addOption('g', "gamma", "weigth to derivatives of images (default: 10)", &o_gamma);
  opt.addOption('s', "scale", "scale [0.5-0.99] for multiscale approach (default: 0.7)", &o_scale);

  opt.parse();          */

class GetOpt {
 protected:
  struct Option {
    enum Type { SWITCH, OPTION, ARGUMENT, OPTIONAL };
    Type type;
    char o;
    QString name;
    QString description;
    QString *value;
    bool *b;
  };
  bool unlimitedArgs;
  QList<Option> options;

 public:
  QString appname;          //application name
  QString help;             //help text
  QStringList args;         //original argument vector
  QStringList arguments;    //arbitrary long list of arguments are store here if unlimitedArgs is true

  GetOpt(): unlimitedArgs(false) {}
  GetOpt(int argc, char *argv[] );
  GetOpt(const QStringList &a);

  //add an option without a value
  void addSwitch(char s, const QString &longname, const QString &description, bool *b );

  //add a valued option (v will be left untouched if the option is not given)
  void addOption(char s, const QString &longname, const QString &description, QString *v );

  //add an argument
  void addArgument(const QString &name, const QString &description, QString *v);

  //add an optional agrument
  void addOptionalArgument(const QString &name, const QString &description, QString *v);
 
  //allow an unlimited number of optional arguments
  void allowUnlimitedArguments(bool allow) {
    unlimitedArgs = allow;
  }
  //set help if someone uses -h or --help option
  void setHelp(QString &_help) { help = _help; }
  
  void parse();

  //return usage string
  QString usage();
  //return argv[0]
  QString &applicationName();

 protected:
  //parses and return true on success
  bool parse(QString &error);
  //return options or switch
  bool findOption(char c, Option &option);
  //return any named argument
  bool findArg(const QString &name, Option &option);
  QString formatDesc(QString desc, int len);
};

#endif
