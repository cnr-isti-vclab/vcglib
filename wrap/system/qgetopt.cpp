
#include <assert.h>
#include "qgetopt.h"

#include <iostream>
using namespace std;

GetOpt::GetOpt(int argc, char *argv[] ) {
  appname = argv[0];
  for(int i = 1; i < argc; i++)
    args.push_back(argv[i]);
}

GetOpt::GetOpt(const QStringList &a): args(a) {
  appname = a[0];
  args = a;
  args.pop_front();
}

 //add an option without a value
void GetOpt::addSwitch(char s, const QString &name, const QString &description, bool *b ) {
  Option option;
  assert(!findOption(s, option));
  assert(!findArg(name, option));
  option.type = Option::SWITCH;
  option.o = s;
  option.name = name;
  option.description = description;
  option.b = b;
  options.push_back(option);
}

  //add a valued option (v will be left untouched if the option is not given)
void GetOpt::addOption(char s, const QString &name, const QString &description, QString *v ) {
  Option option;
  assert(!findOption(s, option));
  assert(!findArg(name, option));
  option.type = Option::OPTION;
  option.o = s;
  option.name = name;
  option.description = description;
  option.value = v;
  options.push_back(option);
}
  //add an argument
void GetOpt::addArgument(const QString &name, const QString &description, QString *v) {
  Option option;
  assert(!findArg(name, option));
  option.type = Option::ARGUMENT;
  option.name = name;
  option.description = description;
  option.value = v;
  options.push_back(option);
}
  //add an optional agrument
void GetOpt::addOptionalArgument(const QString &name, const QString &description, QString *v) {
  Option option;
  assert(!findArg(name, option));
  option.type = Option::OPTIONAL;
  option.name = name;
  option.description = description;
  option.value = v;
  options.push_back(option);
}
  //return application name
QString &GetOpt::applicationName() {
  return appname;
}

  //return usage string
QString GetOpt::usage() {
  QString u = "Usage: " + appname;
  //arguments
  bool has_optionals = false;
  bool has_options = false;
  for(int i = 0; i < options.size(); i++) {
    if(options[i].type == Option::OPTION) has_options = true;
    if(options[i].type == Option::OPTIONAL) has_optionals = true;
    if(options[i].type != Option::ARGUMENT) continue;
    u += " <" + options[i].name + ">";
  }
  //optional arguments
  if(has_optionals) {
    u += " [";
    for(int i = 0; i < options.size(); i++) {
      if(options[i].type != Option::OPTIONAL) continue;
      u += " <" + options[i].name + ">";
    }
    u += "]";
  }
  if(has_options) {
    u += " [-";
    for(int i = 0; i < options.size(); i++) {
      if(options[i].type != Option::OPTION) continue;
      u += options[i].o;
    }
    u += "]";
  }
  u += "\n\n";
  //compute maxlen:
  int maxlen = 0;
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    int len = o.name.size() + 2;
    switch(o.type) {
      case Option::ARGUMENT: 
      case Option::OPTIONAL: break;
      case Option::SWITCH:   len += 5; break;
      case Option::OPTION:   len += 16; break;
      default: break;
    }
    if(len > maxlen) maxlen = len;
  }
  //print options and arguments in the given order
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    QString line = "";
    switch(o.type) {
      case Option::ARGUMENT:
      case Option::OPTIONAL: line += o.name; break;
      case Option::SWITCH:   line += "-" + QString(o.o) + " --" + o.name; break;
      case Option::OPTION:   line += "-" + QString(o.o) + " <val> --" + o.name + "=<val>"; break;
      default: break;
    }
    QString blank = "";
    blank.resize(maxlen - line.size());
    blank.fill(' ');
    line += blank + formatDesc(o.description, maxlen) + "\n"; 
    u += line;
  }
  return u;
}

void GetOpt::parse() {
  QString error;
  if(!parse(error)) {
    cerr << qPrintable(error) << endl << endl << qPrintable(usage()) << endl << endl;
    exit(0);
  }
}

bool GetOpt::parse(QString &error) {
  for(int i = 0; i < args.size(); i++) {
    QString arg = args[i];
    if(args[i] == "-h" || args[i] == "--help") {
      cout << qPrintable(usage()) << endl << qPrintable(help) << endl;
      exit(0);
    }
    //long option
    if(arg.startsWith( QString::fromLatin1( "--" ) ) ) {
      arg = arg.mid( 2 );
      if(arg.isEmpty()) {
        error = "'--' feature not supported, yet";
        return false;
      }
      // split key=value style arguments
      int equal = arg.indexOf( '=' );
      QString val;
      if(equal >= 0) {
        arg = arg.left(equal);
        val = arg.mid(equal + 1);
        if(val.isEmpty()) {
          error = "Emtpy value for option '--" + arg + "'";
          return false;
        }
      }
      if(arg.isEmpty()) {
        error = "Option long name missing: '--=" + val + "' is not a valid option";
        return false;
      }
      Option o;
      if(!findArg(arg, o) || (o.type != Option::OPTION && o.type != Option::SWITCH)) {
        error = "Unknown option: '--" + arg + "'";
        return false;
      }
      *(o.value) = val;

    //option
    } else if( arg[0] == '-' ) {
      if(arg.size() != 2) {
        error = "Invalid option: " + arg;
        return false;
      }
      Option o;
      if(!findOption(arg[1].toAscii(), o)) {
         error = "Unknown option: '" + arg + "'";
        return false;
      }
      if(o.type == Option::SWITCH) {
        *(o.b) = true;
      } else { //OPTION
        i++;
        arg = args[i];
        if(i == args.size() || arg[0] == '-') {
          error = "Missing argument after option '" + arg + "'";
          return false;
        }
        *(o.value) = arg;
      }
    //argument
    } else {
      arguments.push_back(arg);
    } 
  }
  //test arguments
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    if(o.type != Option::ARGUMENT) continue;
    if(arguments.isEmpty()) {
      error = "Too few arguments, could not parse argument '" + o.name + "'";
      return false;
    }
    *(o.value) = arguments.front();
    arguments.pop_front();
  }
  //test arguments
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    if(o.type != Option::ARGUMENT) continue;
    if(arguments.isEmpty()) break;
    *(o.value) = arguments.front();
    arguments.pop_front();
  }
  if(!arguments.isEmpty() && !unlimitedArgs) {
    error = "Too many arguments";
    return false;
  }
  return true;
}

bool GetOpt::findOption(char c, Option &option) {
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    if(o.type != Option::OPTION && o.type != Option::SWITCH) continue;
    if(o.o == c) {
      option = o;
      return true;
    }
  }
  return false;
}

bool GetOpt::findArg(const QString &name, Option &option) {
  for(int i = 0; i < options.size(); i++) {
    Option &o = options[i];
    if(o.name == name) {
      option = o;
      return true;
    }
  }
  return false;
}

QString GetOpt::formatDesc(QString desc, int len) {
  QString output;
  //search for last space before 79 - len characters
  while(!desc.isEmpty()) {
    int pos = desc.lastIndexOf(" ", 79 - len);
    if(pos == -1) {
      output += desc;
      break;
    }
    output += desc.left(pos) + "\n";
    QString blank;
    blank.resize(len);
    blank.fill(' ');
    output += blank;
    desc = desc.mid(pos+1);
  }
  return output;
}