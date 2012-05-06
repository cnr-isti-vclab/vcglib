#ifndef GLW_BOOKKEEPING_H
#define GLW_BOOKKEEPING_H

#include "./common.h"

namespace glw
{

class Context;

namespace detail
{

template <typename TObject>
class SharedObjectBase
{
	friend class Context;

	public:

		typedef void                      BaseType;
		typedef SharedObjectBase<TObject> ThisType;
		typedef TObject                   ObjectType;

		bool isNull(void) const
		{
			return (this->m_object == 0);
		}

		void setNull(void)
		{
			this->m_object = 0;
		}

		const ObjectType * object(void) const
		{
			return this->m_object;
		}

		ObjectType * object(void)
		{
			return this->m_object;
		}

		void ref(void)
		{
			this->m_refCount++;
		}

		void unref(void)
		{
			GLW_ASSERT(this->m_refCount > 0);
			this->m_refCount--;
			if (this->m_refCount == 0)
			{
				if (this->m_object != 0)
				{
					this->signalDestruction();
				}
				delete this;
			}
		}

	protected:

		Context *    m_context;
		ObjectType * m_object;
		int          m_refCount;

		SharedObjectBase(Context * Context, ObjectType * Object)
			: m_context  (Context)
			, m_object   (Object)
			, m_refCount (0)
		{
			GLW_ASSERT(this->m_context != 0);
			GLW_ASSERT(this->m_object  != 0);
		}

	private:

		SharedObjectBase(const ThisType & other);
		ThisType & operator = (const ThisType & other);

		inline void signalDestruction(void);
};

template <typename TObject>
class SharedObject;

template <typename TObject, typename TObjectBase>
class SharedObjectTraits
{
	public:

		typedef void                                     BaseType;
		typedef SharedObjectTraits<TObject, TObjectBase> ThisType;
		typedef TObject                                  ObjectType;
		typedef TObjectBase                              ObjectBaseType;
		typedef SharedObject<ObjectBaseType>             SharedObjectBaseType;
};

template <typename TObject>
class SharedObjectTraits<TObject, void>
{
	public:

		typedef void                                     BaseType;
		typedef SharedObjectTraits<TObject, void>        ThisType;
		typedef TObject                                  ObjectType;
		typedef void                                     ObjectBaseType;
		typedef SharedObjectBase<ObjectType>             SharedObjectBaseType;
};

template <typename TObject>
class SharedObject : public SharedObjectTraits<TObject, typename TObject::BaseType>::SharedObjectBaseType
{
	friend class Context;

	public:

		typedef typename SharedObjectTraits<TObject, typename TObject::BaseType>::SharedObjectBaseType BaseType;
		typedef SharedObject<TObject>                                                                  ThisType;
		typedef TObject                                                                                ObjectType;

		const ObjectType * object(void) const
		{
			return static_cast<const ObjectType *>(BaseType::object());
		}

		ObjectType * object(void)
		{
			return static_cast<ObjectType *>(BaseType::object());
		}

	protected:

		SharedObject(Context * Context, ObjectType * Object)
			: BaseType(Context, Object)
		{
			;
		}

	private:

		SharedObject(const ThisType & other);
		ThisType & operator = (const ThisType & other);
};

template <typename TObject>
class SharedObjectBinding
{
	friend class Context;

	public:

		typedef void                          BaseType;
		typedef SharedObjectBinding<TObject>  ThisType;
		typedef TObject                       ObjectType;
		typedef ObjectType                    UnsafeType;

		~SharedObjectBinding(void)
		{
			this->detach();
		}

		bool isNull(void) const
		{
			if (this->m_shared == 0) return true;
			return this->m_shared->isNull();
		}

		void setNull(void)
		{
			this->m_shared = 0;
		}

		const ObjectType * object(void) const
		{
			GLW_ASSERT(!this->isNull());
			const ObjectType * obj = this->m_shared->object();
			obj->setBinding(this->m_target, this->m_unit);
			return obj;
		}

		ObjectType * object(void)
		{
			GLW_ASSERT(!this->isNull());
			ObjectType * obj = this->m_shared->object();
			obj->setBinding(this->m_target, this->m_unit);
			return obj;
		}

		GLenum target(void) const
		{
			return this->m_target;
		}

		GLint unit(void) const
		{
			return this->m_unit;
		}

		void ref(void)
		{
			this->m_refCount++;
		}

		void unref(void)
		{
			GLW_ASSERT(this->m_refCount > 0);
			this->m_refCount--;
			if (this->m_refCount == 0)
			{
				delete this;
			}
		}

	protected:

		typedef SharedObject<ObjectType> SharedObjectType;

		SharedObjectBinding(SharedObjectType * shared, GLenum target, GLint unit)
			: m_shared   (0)
			, m_refCount (0)
			, m_target   (target)
			, m_unit     (unit)
		{
			this->attach(shared);
		}

	private:

		SharedObjectType * m_shared;
		int                m_refCount;
		GLenum             m_target;
		GLint              m_unit;

		SharedObjectBinding(const ThisType & other);
		ThisType & operator = (const ThisType & other);

		void attach(SharedObjectType * shared)
		{
			this->detach();
			this->m_shared = shared;
			if (this->m_shared != 0)
			{
				this->m_shared->ref();
			}
		}

		void detach(void)
		{
			if (this->m_shared == 0) return;
			this->m_shared->unref();
			this->m_shared = 0;
		}
};

template <typename TObject>
class SafeHandleBase
{
	friend class Context;

	public:

		typedef void                          BaseType;
		typedef SafeHandleBase<TObject>       ThisType;
		typedef TObject                       ObjectType;
		typedef typename ObjectType::SafeType SafeType;

		SafeHandleBase(void)
			: m_shared(0)
		{
			;
		}

		SafeHandleBase(const ThisType & other)
			: m_shared(0)
		{
			this->attach(other.shared());
		}

		~SafeHandleBase(void)
		{
			this->detach();
		}

		bool isNull(void) const
		{
			if (this->m_shared == 0) return true;
			return this->m_shared->isNull();
		}

		void setNull(void)
		{
			this->detach();
		}

		const SafeType * operator -> (void) const
		{
			GLW_ASSERT(!this->isNull());
			return this->m_shared->object();
		}

		SafeType * operator -> (void)
		{
			GLW_ASSERT(!this->isNull());
			return this->m_shared->object();
		}

		ThisType & operator = (const ThisType & other)
		{
			this->attach(other.shared());
			return (*this);
		}

	protected:

		typedef SharedObject<ObjectType> SharedObjectType;

		SafeHandleBase(SharedObjectType * shared)
			: m_shared(0)
		{
			this->attach(shared);
		}

		const ObjectType * object(void) const
		{
			if (this->m_shared == 0) return 0;
			return this->m_shared->object();
		}

		ObjectType * object(void)
		{
			if (this->m_shared == 0) return 0;
			return this->m_shared->object();
		}

		SharedObjectType * shared(void) const
		{
			return this->m_shared;
		}

	private:

		SharedObjectType * m_shared;

		void attach(SharedObjectType * shared)
		{
			this->detach();
			this->m_shared = shared;
			if (this->m_shared != 0)
			{
				this->m_shared->ref();
			}
		}

		void detach(void)
		{
			if (this->m_shared == 0) return;
			this->m_shared->unref();
			this->m_shared = 0;
		}
};

template <typename TObject>
class SafeHandle;

template <typename TObject, typename TObjectBase>
class SafeHandleTraits
{
	public:

		typedef void                                   BaseType;
		typedef SafeHandleTraits<TObject, TObjectBase> ThisType;
		typedef TObject                                ObjectType;
		typedef TObjectBase                            ObjectBaseType;
		typedef SafeHandle<ObjectBaseType>             SafeHandleBaseType;
};

template <typename TObject>
class SafeHandleTraits<TObject, void>
{
	public:

		typedef void                                   BaseType;
		typedef SafeHandleTraits<TObject, void>        ThisType;
		typedef TObject                                ObjectType;
		typedef void                                   ObjectBaseType;
		typedef SafeHandleBase<ObjectType>             SafeHandleBaseType;
};

template <typename TObject>
class SafeHandle : public SafeHandleTraits<TObject, typename TObject::BaseType>::SafeHandleBaseType
{
	friend class Context;

	public:

		typedef typename SafeHandleTraits<TObject, typename TObject::BaseType>::SafeHandleBaseType BaseType;
		typedef SafeHandle<TObject>                                                                ThisType;
		typedef TObject                                                                            ObjectType;
		typedef typename ObjectType::SafeType                                                      SafeType;

		SafeHandle(void)
			: BaseType()
		{
			;
		}

		SafeHandle(const ThisType & other)
			: BaseType(other)
		{
			;
		}

		const SafeType * operator -> (void) const
		{
			return dynamic_cast<const SafeType *>(BaseType:: operator ->());
		}

		SafeType * operator -> (void)
		{
			return dynamic_cast<SafeType *>(BaseType:: operator ->());
		}

		/*
		ThisType & operator = (const ThisType & other)
		{
			this->attach(other.shared());
			return (*this);
		}
		*/

		operator bool (void) const
		{
			return !this->isNull();
		}

	protected:

		typedef SharedObject<ObjectType> SharedObjectType;

		SafeHandle(SharedObjectType * shared)
			: BaseType(shared)
		{
			;
		}

		const ObjectType * object(void) const
		{
			return static_cast<const ObjectType *>(BaseType::object());
		}

		ObjectType * object(void)
		{
			return static_cast<ObjectType *>(BaseType::object());
		}

		SharedObjectType * shared(void) const
		{
			return static_cast<SharedObjectType *>(BaseType::shared());
		}
};

template <typename TObject>
class UnsafeHandle
{
	friend class Context;

	public:

		typedef void                  BaseType;
		typedef UnsafeHandle<TObject> ThisType;
		typedef TObject               ObjectType;
		typedef ObjectType            UnsafeType;

		UnsafeHandle(void)
			: m_shared(0)
		{
			;
		}

		UnsafeHandle(const ThisType & other)
			: m_shared(0)
		{
			this->attach(other.shared());
		}

		~UnsafeHandle(void)
		{
			this->detach();
		}

		bool isNull(void) const
		{
			if (this->m_shared == 0) return true;
			return this->m_shared->isNull();
		}

		void setNull(void)
		{
			this->detach();
		}

		const UnsafeType * operator -> (void) const
		{
			GLW_ASSERT(!this->isNull());
			return this->m_shared->object();
		}

		UnsafeType * operator -> (void)
		{
			GLW_ASSERT(!this->isNull());
			return this->m_shared->object();
		}

		ThisType & operator = (const ThisType & other)
		{
			this->attach(other.shared());
			return (*this);
		}

		operator bool (void) const
		{
			return !this->isNull();
		}

	protected:

		typedef SharedObjectBinding<ObjectType> SharedObjectBindingType;

		UnsafeHandle(SharedObjectBindingType * shared)
			: m_shared(0)
		{
			this->attach(shared);
		}

		const ObjectType * Object(void) const
		{
			if (this->m_shared == 0) return true;
			return this->m_shared->object();
		}

		ObjectType * Object(void)
		{
			if (this->m_shared == 0) return true;
			return this->m_shared->object();
		}

	private:

		SharedObjectBindingType * m_shared;

		void attach(SharedObjectBindingType * shared)
		{
			this->detach();
			this->m_shared = shared;
			if (this->m_shared != 0)
			{
				this->m_shared->ref();
			}
		}

		void detach(void)
		{
			if (this->m_shared == 0) return;
			this->m_shared->unref();
			this->m_shared = 0;
		}

		SharedObjectBindingType * shared(void) const
		{
			return this->m_shared;
		}
};

};

};

#endif // GLW_BOOKKEEPING_H
