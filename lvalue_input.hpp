#ifndef LVALUE_INPUT_HPP__
#define LVALUE_INPUT_HPP__
// http://stackoverflow.com/a/9089311/146041
namespace lvalue_input {
template <typename T>
class in {
	const T*	x;
public:
			in(const volatile T &x_) : x(const_cast<const T*>(&x_)) {}
	const T*	operator-> () const { return x; }
	const T&	operator*  () const { return *x; }
};
} // namespace lvalue_input
#endif
