configure_file (generator_script.m.in ${CMAKE_CURRENT_SOURCE_DIR}/generator_script.m)

foreach (TARGET_Q ${LIST_TARGET_Q})
  foreach (TARGET_N ${LIST_TARGET_N})
    foreach (TARGET_W ${LIST_TARGET_W})
      set (local_filename "test_data_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W}_step0_span0.0.mat")

      configure_file (
        data_generator_qXX_NXX_wXX_step0_0pi_1.m.in
        ${CMAKE_CURRENT_SOURCE_DIR}/data_generator_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W}_step0_0pi_1.m
      )

      # If MATLAB or Octave have been found, generate test vectors. Just create dummy dependencies otherwise.
      if ((OCTAVE_VERSIONS OR MATLAB_RELEASES) AND ENABLE_DATA_AUTOGENERATION)
        add_custom_command (
          OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename}
          COMMAND ${mat_interpreter} ${mat_generator_args} "data_generator_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W}_step0_0pi_1"
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/save_test_vectors.m
                  ${CMAKE_CURRENT_SOURCE_DIR}/compute_ts_score.m
                  ${CMAKE_CURRENT_SOURCE_DIR}/generate_noisy_sequence.m
                  ${CMAKE_CURRENT_SOURCE_DIR}/corr_abs_max.m
                  ${CMAKE_CURRENT_SOURCE_DIR}/pn_sequences.mat
                  ${CMAKE_CURRENT_SOURCE_DIR}/parameters_20210903.mat
          USES_TERMINAL
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          VERBATIM
        )
      endif ()

      add_custom_target (gen_data_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W} DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename})
      unset (${local_filename})
    endforeach ()
  endforeach ()
endforeach ()

set (local_filename "test_data_q64_N60_w4_step4_span0.5.mat")

# If MATLAB or Octave have been found, generate test vectors. Just create dummy dependencies otherwise.
if ((OCTAVE_VERSIONS OR MATLAB_RELEASES) AND ENABLE_DATA_AUTOGENERATION)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename}
    COMMAND
      ${mat_interpreter} ${mat_generator_args}
      "load('parameters_20210903.mat', 'PN64', 'best_N'); fnm = save_test_vectors(PN64, best_N, pi / 2, 4, 4, 0); fprintf('%s generated.\\n', fnm)"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/save_test_vectors.m
            ${CMAKE_CURRENT_SOURCE_DIR}/compute_ts_score.m
            ${CMAKE_CURRENT_SOURCE_DIR}/generate_noisy_sequence.m
            ${CMAKE_CURRENT_SOURCE_DIR}/corr_abs_max.m
            ${CMAKE_CURRENT_SOURCE_DIR}/parameters_20210903.mat
    USES_TERMINAL
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
endif ()

add_custom_target (gen_data_q64_N60_w4 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename})
unset (${local_filename})

set (local_filename "test_data_q64_N60_w34_step8_span2.0.mat")

# If MATLAB or Octave have been found, generate test vectors. Just create dummy dependencies otherwise.
if ((OCTAVE_VERSIONS OR MATLAB_RELEASES) AND ENABLE_DATA_AUTOGENERATION)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename}
    COMMAND
      ${mat_interpreter} ${mat_generator_args}
      "load('parameters_20210903.mat', 'PN64', 'best_N'); fnm = save_test_vectors(PN64, best_N, pi * 2, 34, 8, 0); fprintf('%s generated.\\n', fnm)"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/save_test_vectors.m
            ${CMAKE_CURRENT_SOURCE_DIR}/compute_ts_score.m
            ${CMAKE_CURRENT_SOURCE_DIR}/generate_noisy_sequence.m
            ${CMAKE_CURRENT_SOURCE_DIR}/corr_abs_max.m
            ${CMAKE_CURRENT_SOURCE_DIR}/parameters_20210903.mat
    USES_TERMINAL
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
endif ()

add_custom_target (gen_data_q64_N60_w34 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${local_filename})
unset (${local_filename})
