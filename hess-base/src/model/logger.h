/**************************************************************************
 * This file is part of the Hess project
 * Copyright (C) 2023 Entroforce LLC
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 **************************************************************************/

#pragma once

#include <string>
#include <iostream>
#include <mutex>

#include "exception.h"

namespace hess {
  extern std::mutex Mutex;
  // Singleton class

  class Logger {
  public:

    static Logger& get_instance() {
      std::lock_guard<std::mutex> myLock(Mutex);
      if (!instance) {
        instance = new Logger();
      }
      return *instance;
    }

    static FILE* get_stream() {
      if (log_file == NULL)
        return stderr;
      return log_file;
    }

    static void set_stream(FILE* stream) {
      log_file = stream;
    }

    static void set_path(const char *path) {
      if ((log_file = fopen(path, "w")) == NULL) {
        throw HessException("error opening logfile" + std::string(path));
      }
    }
  private:
    Logger() = default;
    ~Logger() = default;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    static FILE *log_file;
    static Logger* instance;
  };
}
