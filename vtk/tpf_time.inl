#include "tpf_time.h"

#include "../log/tpf_log.h"

#include "vtkDataObject.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        template <typename float_t>
        inline std::vector<float_t> get_timesteps(vtkInformation* info)
        {
            try
            {
                auto time_steps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
                auto time_range = info->Get(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

                if (time_range != NULL)
                {
                    std::vector<float_t> timesteps;
                    std::size_t num_time_steps = 0;

                    for (; time_steps[num_time_steps] < time_range[1]; ++num_time_steps)
                    {
                        timesteps.push_back(static_cast<float_t>(time_steps[num_time_steps]));
                    }

                    timesteps.push_back(static_cast<float_t>(time_steps[++num_time_steps]));

                    return timesteps;
                }
                else
                {
                    return std::vector<float_t>();
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error getting time steps."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error getting time steps."));
            }
        }

        template <typename float_t>
        inline std::pair<float_t, std::size_t> get_timestep(vtkInformation* info, vtkDataSet* alg)
        {
            try
            {
                auto time_steps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
                auto time_range = info->Get(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

                if (time_range != NULL && time_range[0] != time_range[1])
                {
                    float_t current_time_step = static_cast<float_t>(alg->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP()));
                    std::size_t current_time_step_index = 0;

                    float_t min_time_step_diff = std::abs(static_cast<float_t>(time_steps[0]) - current_time_step);
                    std::size_t num_time_steps = 1;

                    for (; time_steps[num_time_steps] < time_range[1]; ++num_time_steps)
                    {
                        if (std::abs(static_cast<float_t>(time_steps[num_time_steps]) - current_time_step) < min_time_step_diff)
                        {
                            min_time_step_diff = std::abs(static_cast<float_t>(time_steps[num_time_steps]) - current_time_step);
                            current_time_step_index = num_time_steps;
                        }
                    }

                    ++num_time_steps;

                    if (std::abs(static_cast<float_t>(time_steps[num_time_steps]) - current_time_step) < min_time_step_diff)
                    {
                        min_time_step_diff = std::abs(static_cast<float_t>(time_steps[num_time_steps]) - current_time_step);
                        current_time_step_index = num_time_steps;
                    }

                    return std::make_pair(static_cast<float_t>(time_steps[current_time_step_index]), static_cast<std::size_t>(current_time_step_index));
                }
                else
                {
                    return std::make_pair(static_cast<float_t>(0.0L), static_cast<std::size_t>(0));
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error getting current time step."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error getting current time step."));
            }
        }

        template <typename float_t>
        inline float_t get_timestep(vtkInformation* info, std::size_t timestep_index)
        {
            try
            {
                auto time_steps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

                return static_cast<float_t>(time_steps[timestep_index]);
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error getting time step."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error getting time step."));
            }
        }

        inline bool is_first_timestep(vtkInformation* info, std::size_t timestep_index)
        {
            try
            {
                auto time_steps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
                auto time_range = info->Get(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

                if (time_range != NULL && time_range[0] != time_range[1])
                {
                    return time_steps[timestep_index] == time_range[0];
                }
                else
                {
                    return true;
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error checking for first time step."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error checking for first time step."));
            }
        }

        inline bool is_last_timestep(vtkInformation* info, std::size_t timestep_index)
        {
            try
            {
                auto time_steps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
                auto time_range = info->Get(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

                if (time_range != NULL && time_range[0] != time_range[1])
                {
                    return time_steps[timestep_index] == time_range[1];
                }
                else
                {
                    return true;
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error checking for last time step."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error checking for last time step."));
            }
        }

        template <typename float_t>
        inline float_t get_timestep_delta(vtkInformation* info, vtkDataSet* alg)
        {
            try
            {
                std::pair<float_t, std::size_t> timestep = get_timestep<float_t>(info, alg);

                if (is_first_timestep(info, timestep.second) && is_last_timestep(info, timestep.second))
                {
                    return static_cast<float_t>(1.0L);
                }

                if (is_last_timestep(info, timestep.second))
                {
                    return timestep.first - get_timestep<float_t>(info, timestep.second - 1);
                }
                else
                {
                    return get_timestep<float_t>(info, timestep.second + 1) - timestep.first;
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error getting time step delta."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error getting time step delta."));
            }
        }
    }
}