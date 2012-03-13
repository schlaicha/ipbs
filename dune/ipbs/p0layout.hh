// layout for codim0 data
#ifndef _P0LAYOUT_H
#define _P0LAYOUT_H
template <int dim>
struct P0Layout
{
  bool contains(Dune::GeometryType gt)
  {
    if ( gt.dim() == dim ) return true;
    return false;
  }
};
#endif
