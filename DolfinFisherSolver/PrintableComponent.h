#ifndef PRINTABLECOMPONENT_H
#define PRINTABLECOMPONENT_H

class PrintableComponent{
	// must be overwritten by component
	// returns string of component details
	virtual std::string asString() = 0;
};

#endif
