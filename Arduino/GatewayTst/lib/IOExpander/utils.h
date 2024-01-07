// IO expander support
#include <PCA9554.h>

extern PCA9554 io;


/*!
    \brief reset the modules connected with IO expander, such as SX1278 and camera
    \param void reset the default pins (pin 6 and pin 7)
    \return no reture
*/
void IOreset(void);