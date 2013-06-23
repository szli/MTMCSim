#ifndef _SPTRWRAPPER_H_INCLUDED
#define _SPTRWRAPPER_H_INCLUDED

#include <boost/shared_ptr.hpp>

template<class T>
	class sptrWrapper
	{
	public:
		/** The default constructor of shared_ptr<T> will be invoked and sptr points to nothing.*/
		sptrWrapper() 
		{}

		sptrWrapper(const T& inner)
		{
			sptr = inner.clone();
		}
		
		/** No action required. When sptr is destructed, it should destruct the T object because it is the only shared_ptr pointing to the T object.*/
		~sptrWrapper()
		{}
		
		sptrWrapper(const sptrWrapper<T>& original)
		{
			if (original.sptr !=0)
			DataPtr = original.DataPtr->clone();
			else
			DataPtr=0;
		}
Wrapper& operator=(const Wrapper<T>& original)
{
if (this != &original)
{
if (DataPtr!=0)
delete DataPtr;
DataPtr = (original.DataPtr !=0) ?
original.DataPtr->clone() : 0;
}
return *this;
}
T& operator*()
{
return *DataPtr;
}
const T& operator*() const
{
return *DataPtr;
}
const T* const operator->() const
{
return DataPtr;
}


	private:
		shard_ptr<T> sptr;


#endif