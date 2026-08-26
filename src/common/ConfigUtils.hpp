/** \file ConfigUtils.hpp
 * \brief Shared hciReduce configuration-loading helpers.
 * \author Jared R. Males
 */

#ifndef ConfigUtils_hpp
#define ConfigUtils_hpp

#include <string>

#include <mx/app/appConfigurator.hpp>
#include <mx/ioutils/stringUtils.hpp>

namespace mx
{
namespace improc
{

/// Load a boolean target supporting bare and attached-value command-line forms.
/** Targets using this loader must be registered with `mx::app::argType::Optional`. A bare command-line occurrence
 * overrides configuration-file values with true, while an attached command-line value is parsed normally.
 */
template <class verboseT>
void loadBoolConfig( mx::app::appConfigurator &config, /**< [in,out] parsed configuration and usage state */
                     bool &value,                      /**< [in,out] default or resolved boolean value */
                     const std::string &name /**< [in] registered target name */ )
{
    const auto found = config.m_targets.find( name );
    if( found == config.m_targets.end() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "boolean configuration target is not registered: " + name );
    }

    mx::app::configTarget &target = found->second;
    target.used = true;
    if( !target.set )
    {
        return;
    }

    if( target.verbosity > 1 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       name + " may be specified only once on the command line" );
    }

    if( target.verbosity > 0 )
    {
        for( std::size_t index = target.sources.size(); index > 0; --index )
        {
            if( target.sources[index - 1] != "command line" )
            {
                continue;
            }
            if( target.values[index - 1] == "true" )
            {
                value = true;
                return;
            }
            if( target.values[index - 1] == "false" )
            {
                value = false;
                return;
            }
            else
            {
                throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                               name + " must be specified as true or false" );
            }
        }

        value = true;
        return;
    }

    if( target.values.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, name + " has no configured boolean value" );
    }
    mx::error_t conversionError = mx::error_t::noerror;
    const bool parsed = mx::ioutils::stoT<bool>( target.values.back(), &conversionError );
    if( conversionError != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, name + " must be specified as true or false" );
    }
    value = parsed;
}

} // namespace improc
} // namespace mx

#endif // ConfigUtils_hpp
