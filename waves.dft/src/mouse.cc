#include "mouse.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

cMouse::cMouse() {
	active = true;
	mouse_fd = 0;
	mouse_ev = new input_event();
	left_button = 0;
	//mouse_st = new mouse_state();
	mouse_fd = open(MOUSE_DEV, O_RDWR | O_NONBLOCK);
	if (mouse_fd > 0) {
		ioctl(mouse_fd, EVIOCGNAME(256), name);
		std::cout << "   Name: " << name << std::endl << std::endl;
		active = true;
		pthread_create(&thread, 0, &cMouse::loop, this);
	}
}

cMouse::~cMouse() {
	if (mouse_fd > 0) {
		active = false;
		pthread_join(thread, 0);
		close(mouse_fd);
	}
	//delete mouse_st;
	delete mouse_ev;
	mouse_fd = 0;
}

void* cMouse::loop(void *obj) {
	while (reinterpret_cast<cMouse *>(obj)->active) reinterpret_cast<cMouse *>(obj)->readEv();
}

void cMouse::readEv() {
	int bytes = read(mouse_fd, mouse_ev, sizeof(*mouse_ev));
	if (bytes > 0) {

	  /*printf("time %ld.%06ld\ttype %d\tcode %d\tvalue %d\n",
		 mouse_ev->time.tv_sec, mouse_ev->time.tv_usec, 
		 mouse_ev->type, mouse_ev->code, mouse_ev->value);*/

		if (mouse_ev->type & EV_REL) {
			mouse_st.axis[mouse_ev->code] = mouse_ev->value;
			//printf("rel\n");
		}
		if (mouse_ev->type & EV_KEY) {
			mouse_st.button[mouse_ev->code] = mouse_ev->value;
			/*printf("key\n");
			printf("Mouse Button [%d] = %d \n",
			BTN_LEFT,mouse_st.button[BTN_LEFT]);*/
			left_button = mouse_ev->value;

		}
	}
}

mouse_state cMouse::getMouseState() {
	mouse_state st = mouse_st;
	mouse_st.axis[0] = 0; mouse_st.axis[1] = 0;
	mouse_st.button[BTN_LEFT]=0;
	return st;
}

bool cMouse::getLeftButton(){
  return left_button;
}
