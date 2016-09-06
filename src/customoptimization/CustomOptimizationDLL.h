#ifdef WIN32
#   ifdef CustomOptimization_EXPORTS
#       define CustomOptimization_API __declspec(dllexport)
#   else
#       define CustomOptimization_API  __declspec(dllimport)
#   endif
#else
#   define CustomOptimization_API
#endif // WIN32