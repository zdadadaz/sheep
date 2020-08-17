#if defined(util)
credit();
#elif defined(DEBIT)
debit();
#else
printerror();
#endif

#include "sheep.c"
bool check_move_valid(int x, int y);
void move(double animal[N][T], int curX, int curY, int t);